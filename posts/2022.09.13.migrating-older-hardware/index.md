@def title = "Migrating to older hardware"
@def date = Date(2022, 09, 13)
@def tags = ["2022", "gemini", "nginx"]

# Migrating to older hardware

_This post is cross-posted from my gemini capsule: `gemini://gemlog.cosroe.com`._

I have been working on migrating the gemini server driving this capsule from a Raspberry Pi 4 to an old Raspberry Pi 2, and transitioning from Titan2 (written in Go) to a server implemented in Zig (and therefore made cooler) in order to flesh out some custom features. The motivation behind this transition is that I need my Raspberry Pi 4 for another project, and with the ongoing electronic component shortage, it is still nearly impossible to acquire new hardware.

The aim is to continue using NGINX as a load balancer and proxy server, however due to the memory and architecture limitations of the Raspberry Pi 2, this can't really be dockerized anymore. I'm also interested in using the NGINX NJS (ECMAScript) stream modules, for which there are unfortunately no armv6 binaries. 

I also have very little experience with Zig, but I am very motivated to learn the language, so this is a good excuse.

## NGINX with NJS from source

To begin with, I installed NGINX with `apt` and used

```bash
sudo nginx -V 2>&1
```

to get all of the build flags used for the armv6 binary. NGINX may subsequently then be removed with `apt` again. 

Next, I cloned the GitHub mirrors of the NGINX source

```bash
git clone https://github.com/nginx/nginx \\
    && git clone https://github.com/nginx/njs \\
    && mv njs nginx/ \\
    && cd nginx
```

This was soon followed by a series of attempts at running the configure script with the build flags pinched from the `apt` binary. As expected, I had a number of missing dependencies, but was able to eventually get everything to work by installing the following

```bash
sudo apt install javascript-common libnginx-mod-stream nodejs libxslt1.1 libxslt1-dev libgeoip-dev
```

I am unsure why I needed to install `libnginx-mod-stream` again, but I received some dynamic linker errors when this wasn't installed. I'm happy to ignore it being there if it's happy to let me run my code c:

The awful configure command that eventually worked:

```bash
./auto/configure \ 
    --with-cc-opt='-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -fPIC -Wdate-time -D_FORTIFY_SOURCE=2' \
    --with-ld-opt='-Wl,-z,relro -Wl,-z,now -fPIC' \
    --prefix=/usr/share/nginx \
    --conf-path=/etc/nginx/nginx.conf \
    --http-log-path=/var/log/nginx/access.log \
    --error-log-path=/var/log/nginx/error.log \
    --lock-path=/var/lock/nginx.lock \
    --pid-path=/run/nginx.pid \
    --modules-path=/usr/lib/nginx/modules \
    --http-client-body-temp-path=/var/lib/nginx/body \
    --http-fastcgi-temp-path=/var/lib/nginx/fastcgi \
    --http-proxy-temp-path=/var/lib/nginx/proxy \
    --http-scgi-temp-path=/var/lib/nginx/scgi \
    --http-uwsgi-temp-path=/var/lib/nginx/uwsgi \
    --with-debug \
    --with-pcre-jit \
    --with-http_ssl_module \
    --with-http_stub_status_module \
    --with-http_realip_module \
    --with-http_auth_request_module \
    --with-http_v2_module \
    --with-http_dav_module \
    --with-http_slice_module \
    --with-threads \
    --with-http_addition_module \
    --with-http_geoip_module=dynamic \
    --with-http_gunzip_module \
    --with-http_gzip_static_module \
    --with-http_sub_module \
    --with-http_xslt_module=dynamic \
    --with-stream=dynamic \
    --with-stream_ssl_module \
    --with-stream_ssl_preread_module \
    --with-mail=dynamic \
    --with-mail_ssl_module \
    --add-dynamic-module=$(pwd)/njs/nginx
```

The output of configure is a Makefile can be built with my whole 1 processors in just under an hour.

```bash
make && sudo make install
```

The `install` command seems to put some things in the wrong place (?), but just symbolically linking the modules directory to wherever it says it needs it seems to be sufficient.

Finally, I link the binary somewhere into my `PATH` and copy the template systemd service file that installing via `apt` seems to leave behind (just `systemctl status --all`), modify the service paths appropriately, and relink it for the daemon. Voila, NGINX avec NJS.

## The zig server

[cozroe source code](https://github.com/fjebaker/cozroe)

The Zig server I am in the process of building (affectionately named cozroe) is based on MasterQ32's `zig-serve` gemini implementation using WolfSSL.

[MasterQ32/zig-serve](https://github.com/MasterQ32/zig-serve)

I wrote a number of small patches before I was able to adequately build it with Zig v0.10, due to some zig shadowing rules that WolfSSL violated. After that, it's pretty trivial to put together a basic file server. My implementation is pretty terrible but it does what it needs to, namely: 

- map `/` to `/index.gmi`.
- sanitize all other paths, and return the file they point to

For my particular raspi, I had a few issues cross-compiling Zig until I found the right CPU architecture -- for posterity, I built with

```bash
zig build -Dtarget=arm-linux-musleabi -Dcpu=arm1176jzf_s -Drelease-small
```

## Monitoring metrics

I am a little curious how many people would actually find this capsule, so the next step was to add a SQLite database to store (hashes of) remote addresses and the time of connection, so I can do a little bit of metric monitoring. With nektro's amazing zigmod, this is only a few zig dependencies away, and soon I was connected to an SQLite database with command line arguments parsing.

[nektro/zigmod](https://github.com/nektro/zigmod)
[Hejsil/zig-clap](https://github.com/Hejsil/zig-clap)
[vrischmann/zig-sqlite](https://github.com/vrischmann/zig-sqlite)

But oh no; NGINX proxy pass for TCP streams does't have facilities for forwarding remote addresses, and so all of my connections are localhosts D:

This is where I tried a number of things to get the remote address to the backend before the backend has finished handling the request:

- use Telegraf to inject straight into the SQLite database. This turned out to be far too heavy for my little raspi and also it turns out logs are written after NGINX has finished processing the full request, which would mean the zig server would be waving goodbye to the response before it was told who sent the initial request.

- pass the remote address as a query parameter in the URL, but this doesn't work because the TLS handshake takes place at the backend and not in NGINX, so I can't modify the URL.

- use the proxy protocol to wrap additional information. Again, a bit tricky without heavily modifying zig-serve, which I don't feel confident enough to do with my current zig understanding.

So instead, this is where I need the ECMAScript module for NGINX, to inject a little bit of custom behaviour when NGINX picks up a stream:

```js
var fs = require('fs');

function validate(s) {
    s.on('upload', async function (data, flags) {
        var time_millis = (new Date()).getTime();
        var info = `${time_millis}|${s.variables.time_local}|${s.remoteAddress}\n`;
        // write to a file
        fs.appendFileSync("/tmp/cozroe.gemini.addr", info);

        s.done(0);
    });
}

export default {validate}
```

I don't know why I called this function validate -- it made sense at the time, but what it does is it gets a stream packet, attaches a callback when the packet has new data for the backend, and then logs two time formats and the remote address to a temporary file, before clearing the packet to continue.

I add this to my configuration file using the `js_access` directive:

```bash
stream {
    # ... 

    js_access sql_logger.validate;

    # gemini server
    server {
        access_log /var/log/nginx/gemini_access.log gemini_log;
        error_log /var/log/nginx/gemini_error.log warn;

        # ...
    }
    
    # ...
}
```

The rest of the configuration is extremely similar to my previous post about NGINX setups for gemini. Next I just have to get cozroe to read it and dump it into the database, which is pretty straightforward.

Currently, this would all fall apart under heavy traffic, but fortunately an easy extension to add is to compare the timestamps in the file with the time read by zig to work out roughly which request is currently being dealt with.

But that's for a much, *much* later period in the future.

And with that, our capsule is now running on older hardware, with more features, and able to count the total visits and unique visits to our little capsule c:
