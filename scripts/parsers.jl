using Dates 
using Base

# global consts
const POSTS_DIRECTORY = "posts"
const POSTS_MAIN = "index.md"

struct BlogPost
    path::String
    author::String
    title::String
    summary::String
    date::DateTime
    tags::Vector{String}
end

Base.isless(p1::BlogPost, p2::BlogPost) = p1.date < p2.date

function BlogPost(f)
    author = pagevar(f, :author)
    title = pagevar(f, :title)
    summary = pagevar(f, :summary)
    date = pagevar(f, :date)


    tags = pagevar(f, :tags)
    if isnothing(tags) || isempty(tags)
        @warn "$f has no tags field."
        tags = String[]
    end

    if isnothing(author)
        error("Post $f has no author")
    end
    if isnothing(title)
        error("Post $f has no title")
    end
    if isnothing(summary)
        @warn("Post $f has no summary")
        summary = "No summary"
    end
    if isnothing(date)
        error("Post $f has no date")
    end

    BlogPost(
        # trim off the `.md` extension
        replace(f, POSTS_MAIN[1:end-3] => ""), author, title, summary, date, tags
    )
end

function _get_posts()
    # get all subdirectories
    subdirs = filter(isdir, readdir(POSTS_DIRECTORY, join=true))
    # filter those which contain a `POSTS_MAIN` file
    post_subdirs = filter(subdirs) do d
        files = readdir(d)
        POSTS_MAIN in files
    end
    fs = joinpath.(post_subdirs, POSTS_MAIN)
    # trim `.md` extension
    fs = [i[1:end-3] for i in fs]
    map(BlogPost, fs)
end

_get_posts_from_files(fs) = map(BlogPost, fs)

function format_summary(b::BlogPost)
    date = Dates.format(b.date, "d u Y")
    selected_tags = format_tag.(b.tags)
    last_index = min(length(selected_tags), 5)
    selected_tags = selected_tags[1:last_index]
    tags = join(selected_tags, ",")
    """
    <div class="post-card">
        <div>
            <a id="post-card-selection" href="/$(b.path)">
                <h3>> $(b.title)</h3> 
                <small style="color: grey;"><i> $(date)</i></small>
                <p style="text-align: right;"> <small id="view-text"> ☞ View post </small> </p>
            </a>
        </div>
    </div>
    """
end

function hfun_posts_chronological(n)
    posts = _get_posts()
    sort!(posts; rev=true)

    last_index = min(parse(Int, first(n)), length(posts))
    posts = posts[1:last_index]

    return join(format_summary.(posts), "\n")
end