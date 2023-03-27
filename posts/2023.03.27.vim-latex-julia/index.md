@def title = "Vim + LaTeX + Julia = <3"
@def date = Date(2023, 03, 27)
@def tags = ["2023", "vim", "julia", "latex"]

# {{fill title}}

I want to make this year the year I finally fully move to using _solely_ vim. My reluctance has been to some extent driven by the ease of tooling in VSCode. In particular the language servers and language extensions are sometimes _too_ good, in that it develops an inability to go without.

I have been using the vim extension in VSCode for a few months. I feel comfortable enough navigating using normal mode, and I have for some time been pushing up against the limits of what the vim extension can do in VSCode. So over the weekend I decided to properly invest time making vim do what I have become so used to in VSCode.

\tableofcontents

This post is to document my own choices, to describe the flavour of vim I am trying to brew, and to hopefully convince friends that the potion is worth drinking.

## Note taking

I try to use the [Zettelkasten](https://en.wikipedia.org/wiki/Zettelkasten) method for note taking. For the past year been using [Dendron](https://www.dendron.so/) for VSCode as a means of note taking. There are a number gripes I have with Dendron however, mostly related to some occasionally buggy formatting in titles (especially when linking titles with LaTeX), no equation / figure references,  and lack of custom environments. 

I'd like to note I am fully aware that some (likely all) of these things could be fixed with some clever customization (or by reading the documentation fully), but that is _clearly beyond my capability_. 

But there are also things that I think Dendron does really well. The flat directory structure, and hierarchical file naming conventions that mirror object orientation is inspired, along with how it automatically manages links and backlinks. The note lookup is likewise extremely good, and is the first feature I want to add to vim.

The idea is to use LaTeX to take notes, using some external program to eventually manage links in the same way that Dendron does, but that is lower priority than being able to quickly lookup and find notes.

### Everyone should be using fzf

A [command line fuzzy finder](https://github.com/junegunn/fzf), fzf is absolutely incredible. I have known about the tool for a while, mainly through [rip grep all](https://github.com/phiresky/ripgrep-all), but only now have I properly started looking into it. The linked repository links to a YouTube video of [@samoshkin](https://github.com/samoshkin) showcasing the features and I would highly [recommend watching it](https://www.youtube.com/watch?v=qgG5Jhi_Els).

Better yet: it has full vim support. That means it brings up a beautiful preview dialog for whatever you're fuzzy finding, and lets you sink the output to whatever command you like (edit it in a new tab, open it in a PDF viewer, delete it, etc).

I added a couple of shortcuts for `:Files`, `:BLines`, etc., and already I am able to start navigating through the file system quicker than I ever could with netrw.

So to mimic the lookup command from Dendron, I added two functions:

```c
function! SearchNoteNames()
    call fzf#run(fzf#wrap({'source': 'ls -1 *.tex', 'sink': 'e', 'options': '--reverse'}))
endfunction
function! SearchNotePDFs()
    call fzf#run(fzf#wrap({'source': 'ls -1 build/*.pdf', 'sink': '!open -a skim', 'options': '--reverse'}))
endfunction
```

The first searches for `.tex` files in the current directory, and lists them. Whichever is selected is then opened in a buffer. The second searches in the local `build` directory for any `.pdf` files, and pipes the selected to be opened with [Skim](https://skim-app.sourceforge.io/). Perfect.

### Vim-TeX and multifile projects

### Type less with UltiSnips

## Julia development

### You complete me COC

### Slipping into Vim-Slime

## Putting it all together

## Coda

So. Farewell then VSCode. You were the text editor Microsoft forked, and that has been sending my "diagnostic and usage" information to far off server ever since. VSCodium was a pall, but vim will forever be the land of milk and honey, and fzf the guide to help us find our way there.
