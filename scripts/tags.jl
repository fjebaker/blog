using Colors
using Random
Random.seed!(42)

@memoize _tag_color(tag) = "#" * hex(HSV(floor(Int, rand() * 1000), 1, 1)) 

function _get_all_tags()
    tags = [i[1] => length(i[2]) for i in globvar("fd_tag_pages")]
    sort!(tags; by=last, rev=true)
    tags
end

function format_tag(tag)
    color = _tag_color(tag)
    """<a href="/tag/$tag"> $tag </a><div class="tag-color-box" style="background-color: $color;"></div>"""
end

@delay function hfun_all_tag_list()
    tags = _get_all_tags()
    fmt_tags = [(format_tag(t), "$count post" * (count == 1 ? "" : "s")) for (t, count) in tags]
    html_tags = ["<li> $(t) ($(count)) </li>" for (t, count) in fmt_tags]
    return "<ul>" * join(html_tags, "\n") * "</ul>"
end

function hfun_format_tag()
    tag = locvar(:fd_tag)
    return format_tag(tag)
end

# needed when building `tags.html`
function hfun_postcard_taglist()
    tag = locvar(:fd_tag)
    rpaths = globvar("fd_tag_pages")[tag]
    posts = _get_posts_from_files(rpaths)
    sort!(posts; rev=true)

    return join(
        format_summary.(posts), "\n"
    )
end

@delay function hfun_maketitle()
    title = locvar(:title)
    author = locvar(:author)
    affiliation = locvar(:affiliation)
    date = locvar(:date)
    tags = locvar(:tags)
    fmt_tags = if !isempty(tags)
        join(format_tag.(tags), ", ")
    else
        ""
    end
    
    return """
    <div class="post-title">
    <h1> $title </h1>
    <small style="margin-right: 1em;"> <i> $(fmt_tags) </i> </small>
    <div style="margin-top: 0; width : 100%;">
        <div style="float: left; text-align: left;"> 
            <small> <i> Author: $author ($affiliation)</i> </small> 
        </div>
        <div style="text-align: right;"> 
            <small> <i> $date </i> </small> 
        </div>
    </div>
    </div>
    """
end