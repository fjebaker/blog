<!--
Add here global page variables to use throughout your website.
-->
+++
author = "Fergus Baker"
mintoclevel = 2

# Add here files or directories that should be ignored by Franklin, otherwise
# these files might be copied and, if markdown, processed by Franklin which
# you might not want. Indicate directories by ending the name with a `/`.
# Base files such as LICENSE.md and README.md are ignored by default.
ignore = ["node_modules/"]

# RSS (the website_{title, descr, url} must be defined to get RSS)
generate_rss = false
website_title = "/var/log/user.log"
website_descr = "..."
+++

<!--
Add here global latex commands to use throughout your pages.
-->
\newcommand{\tensor}[3]{#1^{#2}_{\phantom{#2}#3}}
\newcommand{\d}{\text{d}}
\newcommand{\deriv}[2]{\frac{\d#1}{\d#2}}
\newcommand{\pderiv}[2]{\frac{\partial#1}{\partial#2}}
