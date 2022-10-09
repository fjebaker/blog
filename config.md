+++
author = "Fergus Baker"
mintoclevel = 2
prepath = "blog"

ignore = ["node_modules/"]

generate_rss = false
website_title = "/var/log/user.log"
website_descr = "..."
+++

<!--
global latex commands
-->
\newcommand{\tensor}[3]{#1^{#2}_{\phantom{#2}#3}}
\newcommand{\btensor}[3]{\tensor{\mathbf{#1}}{#2}{#3}}
\newcommand{\ddelta}[2]{\delta^{#1}_{#2}}
\newcommand{\frame}[2]{ \frac{\partial}{\partial \tensor{#1}{#2}{} }}
\newcommand{\basis}[2]{ \d \tensor{#1}{#2}{} }
\newcommand{\const}{\text{const}}
\newcommand{\d}{\text{d}}
\newcommand{\e}{\text{e}}
\newcommand{\deriv}[2]{\frac{\d#1}{\d#2}}
\newcommand{\pderiv}[2]{\frac{\partial#1}{\partial#2}}
\newcommand{\sign}{\text{sign}\,}
\newcommand{\iprod}[2]{\langle #1 , #2 \rangle}