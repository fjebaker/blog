using Dates

function hfun_bar(vname)
  val = Meta.parse(vname[1])
  return round(sqrt(val), digits=2)
end

function hfun_m1fill(vname)
  var = vname[1]
  return pagevar("index", var)
end

function _build_summary(page)
    title = pagevar(page, :title)
    description = pagevar(page, :description)
    image = pagevar(page, :image)
    date = Dates.format(
        Date(pagevar(page, :date), dateformat"d/m/y"),
        "d u Y"
    )
    tags = join(("""<a href="tag/$i">$i</a>""" for i in pagevar(page, :tags)), " - ")

    """<div class="row"><div class="container">
    <img class="summary-image" src="$image">
    <a href="$page"><h2>$title</h2></a>
    <p>$description</p>
    <p style="font-size:14px;">
    $date - $tags
    </p>
    <div style="clear: both"></div>
    </div>
    </div>
    """
end

function hfun_post_summaries()
    pages = [
        joinpath("pages", i)[1:end-3] for i in filter(f -> endswith(f, ".md"), readdir("pages"))
    ]
    perm = sortperm([pagevar(p, :date) for p in pages])
    io = IOBuffer()
    write(io, _build_summary(pages[1]))
    if length(perm) > 1
        for i in perm[2:end]
            summary = _build_summary(pages[i])
            write(io, "<hr/>\n" * summary)
        end
    end
    String(take!(io))
end

function lx_baz(com, _)
  # keep this first line
  brace_content = Franklin.content(com.braces[1]) # input string
  # do whatever you want here
  return uppercase(brace_content)
end
