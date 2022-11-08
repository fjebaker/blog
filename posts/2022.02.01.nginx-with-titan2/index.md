
@def title = "Using NGINX with Titan2 to launch capsules"
@def date = Date(2022, 02, 01)
@def tags = ["2022", "gemini", "nginx"]

# Using NGINX with Titan2 to launch capsules

_This post is cross-posted from my gemini capsule: `gemini://gemlog.cosroe.com`._

This is a short guide to self-hosting a Gemini capsule on a Raspberry Pi (optionally with docker). An NGINX image is used to direct traffic to an upstream Titan2 Gemini file server.

[lostleonardo's Titan2 Gemini server](https://gitlab.com/lostleonardo/titan2)

The Raspberry Pi I am planning to use I am already (ab)using with multiple services. Unfortunately, at time of writing there is an electronic components shortage, which has made it extremely difficult (and expensive) to procure more hardware to run these services on -- my motivation therefore lies in creating something I can easily move to a different machine without having to spend time reconfiguring the setup, once/if the shortage is addressed. Consequently, I will be using docker to spawn the NGINX server, but this is not at all required -- the mechanism by which the service is spawned is irrelevant, though for posterity I have included docker scripts at the end for those similarly interested. 

The motivation for using NGINX relates to load-balancing, centralised logging, and to co-exist with HTTP servers. I won't cover how to set those things up in this post.

I am using:
- Raspberry Pi 4 (4GB)
- 64-bit Raspberry Pi OS 

## NGINX setup

The NGINX server is configured through a single nginx.conf file, inspired by panda-roux's NGINX setup

[panda-roux's capsule](gemini://gemini.panda-roux.dev/log/entry)

```
# conf/nginx.conf
events {
    worker_connections 1024;
}
stream {
    # Gemini server

    map $ssl_preread_server_name  $name {
        gemlog.cosroe.com titan2_server;
    }
    
    server {
        listen 1965;
        ssl_preread on;
        proxy_buffer_size 16k;

        # pass to the backend
        proxy_pass $name;
    }

    upstream titan2_server {
        # local IP address
        server $titan_ip_addr:1961;
    }
}
```

Note, "titan_ip_addr" in the "upstream" directive is localhost if the Titan2 instance is running on the same machine (for docker, it will be the local ip address of that machine). It also redirects to port 1961, since NGINX will be listening on 1965.

## Titan2

Titan2 is implemented in go, and so requires a go runtime. Installing this on a Raspberry Pi is straight forward; select the correct binaries, and then fetch:
```
wget https://dl.google.com/go/go1.17.6.linux-arm64.tar.gz -O go.tar.gz
sudo tar -C /usr/local -xzf go.tar.gz

# add env variables
echo "export GOPATH=$HOME/go" >> .bashrc
echo "export PATH=/usr/local/go/bin:$PATH:$GOPATH/bin" >> .bashrc
```

Building the project is then very easy 
```
git clone "https://gitlab.com/lostleonardo/titan2"
cd titan2
go install
```
In order to launch the capsule, we require self-signed certificates for our domain. These can be generated with openssl:
```
mkdir -p certs
openssl req -x509 -newkey rsa:4096 \
    -keyout certs/key.rsa \
    -out certs/cert.pem \
    -days 3650 -nodes \
    -subj "/CN=gemlog.cosroe.com"
```

Finally, we require some content to serve
```
mkdir -p capsule \
    && echo -e "# Hello Gemini\nCapsule launched $(date)" \
    > capsule/index.gmi
```

We are then ready to launch!
```
titan2 -hostname 0.0.0.0 \
    -dir capsule \
    -crt certs/cert.pem \
    -key certs/key.rsa \
    -port 1961
```

## Docker and systemctl

NGINX can be configured as a systemctl service with a simple service file. When using docker (as I am), this may look like
```
# nginx.service
[Unit]
Description=NGINX Service
Requires=docker.service
After=docker.service

[Service]
Type=oneshot
RemainAfterExit=yes
ExecStart=/usr/bin/docker run --name nginx --rm -d \
	-v /home/pi/conf/nginx.conf:/etc/nginx/nginx.conf \
	-p "1965:1965" \
	nginx
ExecStop=/usr/bin/docker stop nginx
TimeoutStartSec=0

[Install]
WantedBy=multi-user.target
```

To ensure Titan2 launches on startup, we can create a simple service file:
```
# titan2.service
[Unit]
Description=titan2 gemini server
After=network.target

[Service]
User=pi
Type=simple
ExecStart=/home/pi/go/bin/titan2 \
    -hostname 0.0.0.0 \
     -dir /home/pi/capsule \
     -crt /home/pi/certs/cert.pem \
     -key /home/pi/certs/key.rsa \
     -port 1961

[Install]
WantedBy=default.target
```

Then we link and enable the services:
```
sudo systemctl link $(pwd)/nginx.service
sudo systemctl link $(pwd)/titan2.service

sudo systemctl enable nginx titan2
```

