%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LaTeX document class: cseethesis.cls
%% Licentiate or doctoral thesis style file,
%% Luleå Univ. of Tech., Sweden.
%%
%% Created by: Johan E. Carlson, November 5, 2002, originally named 
%%             eisthesis.cls
%%
%% Derived from the Standard LaTeX Book document class.
%%
%% This document class is free to use, copy, distribute,
%% and modify, as long as reference to the original author
%% is given.
%%
%% Last modified by: Johan E. Carlson, June 7, 2009.
%% Additional modification by Tosin Adewumi, Oct 9, 2021.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Revisions (see example document for more information)
%% -------------------------------------------------------
%%   April 7, 2004:
%%      Removed '\ifnum' error
%%      Removed conflict between this template and the
%%        amsmath package.
%%   November 17, 2004:
%%      Corrected float-to-text separation error.
%%   November 19, 2004:
%%      Included support for EU logo and pdfLaTeX/LaTeX
%%      class options.
%%   February 22, 2005:
%%      Major bug fixes (see example document)
%%   January 17, 2007:
%%      See example document.
%%   January 19, 2007:
%%      See example document.
%%   January 25, 2009:
%%      See example document
%%   February 2, 2009
%%   March 5, 2009
%%   June 7, 2009
%%   September 1, 2010
%%  Oct 9, 2021
%%  Addition of makepapertobesubmitted which is still in manuscript stage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ******************************************************
% Identification
% ******************************************************
\def\classversion{3.1}
\def\classdate{October 9, 2021}

\typeout{Document class version: \classversion}
\typeout{Last modified: \classdate}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cseethesis}

% ******************************************************
% Declare Options
% ******************************************************
% Do not know how to handle other options and therefore
% just pass them on to 'book'
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}

% ******************************************************
% Process Options
% ******************************************************
\ProcessOptions*

% ******************************************************
% Load Used Packages and Classes
% ******************************************************
\LoadClass{book}
\RequirePackage{graphicx}
\RequirePackage{fancyhdr}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{ifthen}
\RequirePackage[sectionbib]{bibunits}
\usepackage{multicol}
\usepackage{nomencl}
\usepackage[utf8]{inputenc}
\usepackage{float}
\usepackage[acronym]{glossaries}
\usepackage{acronym}
\usepackage{comment}
\usepackage{hyperref}
\usepackage{lipsum}
\usepackage{tikz}
\usepackage{tikz-3dplot}
\usepackage{xcolor}
\usetikzlibrary{arrows.meta, calc, positioning}
%\RequirePackage[latin1]{inputenc}
\RequirePackage{titletoc}

% ******************************************************
% Start defining and redefining old book style stuff
% ******************************************************

\setlength{\topmargin}{0in}
\setlength{\textwidth}{6.25in}
\setlength{\textheight}{8.75in}
\setlength{\evensidemargin}{0mm}%-4.5mm originally % otherwise every other page as messed up margins....... who made this fucking template
\setlength{\oddsidemargin}{0mm}% 5.5mm originally
\setcounter{tocdepth}{1}
\setcounter{topnumber}{2}
\setcounter{totalnumber}{3}




\renewcommand{\thechapter}{\arabic{chapter}}
\renewcommand{\thesection}{\arabic{chapter}.\arabic{section}}
\renewcommand{\thesubsection}{\arabic{chapter}.\arabic{section}.
                              \arabic{subsection}}
\renewcommand{\thesubsubsection}{\arabic{chapter}.\arabic{section}.
                                 \arabic{subsection}.\arabic{subsubsection}}
\renewcommand{\theequation}{\arabic{chapter}.\arabic{equation}}
\renewcommand{\thefigure}{\arabic{chapter}.\arabic{figure}}
\renewcommand{\thetable}{\arabic{chapter}.\arabic{table}}
\renewcommand{\topfraction}{.95}
\renewcommand{\textfraction}{.05}
\renewcommand{\floatpagefraction}{.75}
\renewcommand{\jot}{.6\baselineskip}
\setlength{\intextsep}{24pt}
\setlength{\headheight}{15pt}

\newcommand{\ischapter}{1} % 1=chapters, 0=papers (used in contents page generation)

%----------------------------------------------
% Make chapter headings look nice
%----------------------------------------------
\def\@chapter[#1]#2{
	 \refstepcounter{chapter}
    {\ifthenelse{\equal{\ischapter}{1}}
       {
	          \addcontentsline{toc}{chapter}{{\sc\chaptername}\space%
                \textnormal{\thechapter\space--} \sc #1}%
       }
       {
           \addcontentsline{toc}{chapter}{{\sc\chaptername}\space%
                \textnormal{\thechapter} \sc #1}%
       }
    }
	 \chaptermark{#1}
    \@makechapterhead{#2}%
    \@afterheading%                % No indendation of first row after heading
}

\def\@makechapterhead#1{%
  {\parindent \z@ \raggedleft
    \reset@font
    \vspace*{2cm}
    \hrule
    \vspace*{8\p@}%
    \par
    \Large {\Huge{\sc\@chapapp{}}\textnormal{\space\thechapter}}
    \par\nobreak
    \vspace*{8\p@}%
    \par
    \vspace*{1\p@}%
    \hrule
    \vspace*{8\p@}
    \parbox{12.5cm}
     {
       \begin{flushright}
          \Huge \textnormal{#1}\par\nobreak
       \end{flushright}
     }
    \vskip 50\p@
  }
}

% ------------------------------------------
% Make figure captions \it and \small
% ------------------------------------------
\long\def\@makecaption#1#2{
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\small{\it{#1: #2}}}
  \ifdim \wd\@tempboxa >\hsize
    \small{\it{#1: #2\par}}
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}
  \fi
  \vskip\belowcaptionskip}

% ------------------------------------------
% The \startpreamble command.
% Generates all the stuff before the first
% "real" chapter starts, like Abstract,
% Title page, TOC, etc..
% ------------------------------------------
\newcommand{\startpreamble}[9]
{
  \pagestyle{plain}
  \thispagestyle{empty}

  %--- Make title page
  \setcounter{page}{1}
  \pagenumbering{roman}
  \begin{center}
     \vspace{1cm} \mbox{}\hrulefill \mbox{}\\[2.5cm]
     {\Huge \textbf{#1\\[3cm]}}
     {\LARGE \textbf{#2}}\\[20mm]
     {\large #3} \vfill
     {\large \mbox{}\hrulefill\mbox{}\\[10mm] 
      \textbf{\supervisorstring}\\[3mm]#4\\[12mm]}
     #8
  \end{center}

  \clearpage

  %--- Make page with the \dedication string

  \thispagestyle{plain}
  \vspace*{6cm} \hfill{\emph{#5}}

   %--- Make abstract page
   \chapter*{\hrule\medskip\hfill \sc Abstract\medskip\hrule}
   \noindent #6

   \clearpage

   %--- Include preface

 \setcounter{secnumdepth}{-2}
 \chapter*{\hrule\medskip\hfill \sc Acknowledgments \medskip\hrule}
 \setcounter{secnumdepth}{2}
 \noindent #7
 \clearpage

%--- Make list of figures
{
\let\clearpage\relax
\let\newpage\relax
\let\pagebreak\relax

 \chapter*{\hrule\medskip\hfill \sc List of Figures \medskip\hrule}
 \renewcommand{\listfigurename}{\vspace{-104pt}}
 \begin{multicols}{2}
    \listoffigures
\end{multicols}
} \clearpage

{
\let\clearpage\relax
\let\newpage\relax
\let\pagebreak\relax

%--- Make list of tables
\chapter*{\hrule\medskip\hfill \sc List of Tables \medskip\hrule}
\renewcommand{\listtablename}{\vspace{-104pt}}
\begin{multicols}{2}
    \listoftables
 \end{multicols}
} \clearpage

{
\let\clearpage\relax
\let\newpage\relax
\let\pagebreak\relax

%--- Make list of abbreviations
\chapter*{\hrule\medskip\hfill \sc List of Abbreviations \medskip\hrule}
\begin{multicols}{2}
    \input{acro_list}
\end{multicols}
} \clearpage

{
\let\clearpage\relax
\let\newpage\relax
\let\pagebreak\relax

%--- Make list of symbols
\chapter*{\hrule\medskip\hfill \sc List of Symbols \medskip\hrule}
\renewcommand{\nomname}{\vspace{-104pt}}
\setlength{\nomlabelwidth}{1cm}
\begin{multicols}{2}
    \printnomenclature
\end{multicols}
} \clearpage

   %--- Make Table of Conents
  \renewcommand{\contentsname}{\hrule\vspace*{8pt}\Huge{\hfill{\sc
                Table of Contents}\vspace*{8pt}\hrule\vspace*{-1cm} }}
    
  \titlecontents{chapter}[0.7cm]
    {\vspace{0.3cm}}
    {\contentslabel{1.7cm}}
    {\hspace{-0.7cm}}
    {\hfill\contentspage}
 
   \tableofcontents
   \clearpage
  
    \pagenumbering{arabic}
}

% ------------------------------------------
% Make \startchapters command
% ------------------------------------------
\newcommand{\startchapters}
{
 \clearpage
 \renewcommand{\chaptername}{Chapter}
 \renewcommand{\ischapter}{1}
 \pagestyle{fancy}
 \fancyhf{}
 \fancyhead[RE]{}
 \fancyhead[LO]{}
 \fancyhead[LE]{\thepage} 
 \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% Make \startpapers command
% ------------------------------------------
\newcommand{\startpapers}
 {
  % ------------------------------------------
  % Adjust section headings
  % ------------------------------------------
  \renewcommand\section{\@startsection {section}{1}{-1pt}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\Large\bfseries}}

  \renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\large\bfseries}}
  \renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalsize\bfseries}}

  %\setcounter{secnumdepth}{3}
  \renewcommand{\thesection}{\arabic{section}}
  \renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}
  \renewcommand{\thesubsubsection}{\arabic{section}.\arabic{subsection}.
                                   \arabic{subsubsection}}
  \renewcommand{\theequation}{\arabic{equation}}
  \renewcommand{\thefigure}{\arabic{figure}}
  \renewcommand{\thetable}{\arabic{table}}
  \renewcommand{\chaptername}{Paper}
  \renewcommand{\thechapter}{\textnormal{\Alph{chapter}}}
  \setcounter{chapter}{0}
  \def\bibname{References}
  \renewcommand{\ischapter}{0}
  \pagestyle{fancy}
  \fancyhf{}
  \fancyhead[RE]{}
  \fancyhead[LO]{}
  \fancyhead[LE]{\thepage} 
  \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% Make new part page
% ------------------------------------------
\newcommand{\makepartpage}[1]{%
 \clearpage
 \thispagestyle{plain}
 \vspace*{7cm}
 {\Huge\hfill #1}
 \addcontentsline{toc}{part}{#1}
 \pagestyle{fancy}
 \fancyhf{}
 \fancyhead[RE]{}
 \fancyhead[LO]{}
 \fancyhead[LE]{\thepage} 
 \fancyhead[RO]{\thepage}
 \clearpage
}

% ------------------------------------------
% \makechapter command
% ------------------------------------------
\newcommand{\makechapter}[4][]
{
   % First argument is an abbreviated version
   % of the chapter name, used as page header
   % Second argument is the chapter name as it occurs in the TOC. 
   % The third argument is what is actually printed on the first page
   % of the chapter.

  % ------------------------------------------
  % Adjust section headings
  % ------------------------------------------
  \renewcommand\section{\@startsection {section}{1}{-1pt}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\Large\bfseries}}

  \renewcommand\subsection{\@startsection{subsection}{2}{-6pt}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\large\bfseries}}
  \renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalsize\bfseries}}


   \renewcommand{\thechapter}{\arabic{chapter}}
   \renewcommand{\thesection}{\arabic{chapter}.\arabic{section}} 
   \renewcommand{\thesubsection}{
                 \arabic{chapter}.\arabic{section}.\arabic{subsection}}
   \renewcommand{\theequation}{\arabic{chapter}.\arabic{equation}}
   \renewcommand{\thefigure}{\arabic{chapter}.\arabic{figure}}
   \renewcommand{\thetable}{\arabic{chapter}.\arabic{table}}

   \clearpage
   \thispagestyle{plain}
   \chapter[#3]{#4}
   
   \vspace*{-2\baselineskip}
   \hfill\parbox[t]{.6\textwidth}{%
      \begin{flushright}
		   \textit{#1}
      \end{flushright}
   }

   \fancyhf{}
   \fancyhead[RE]{\sc{\nouppercase{#2}}}
   \fancyhead[LO]{\sc{\nouppercase{\rightmark}}}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% \makepaper command
% ------------------------------------------
\newcommand{\makepaper}[5]
{
   \fancyhead[RE]{}
   \fancyhead[LO]{}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
   \clearpage
   \thispagestyle{plain}
   \chapter[\space]{#2}
   \vspace{25mm}
   \textbf{Authors:}\\ #3
   \vspace{10mm} \\ \noindent \textbf{Reformatted version of paper
    originally published in:} \\ #4
   \vspace{15mm} \\ \noindent{\copyright} #5
   \clearpage
   \fancyhf{}
   \fancyhead[RE]{\sc{\nouppercase{#1}}}
   \fancyhead[LO]{\sc{\nouppercase{\rightmark}}}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% \makepaperaccepted command, same as above
% for submitted, not yet accepted paper
% ------------------------------------------
\newcommand{\makepaperaccepted}[5]
{
   \fancyhead[RE]{}
   \fancyhead[LO]{}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
   \clearpage
   \thispagestyle{plain}
   \chapter[\space]{#2}
   \vspace{25mm}
   \textbf{Authors:}\\ #3
   \vspace{10mm} \\ \noindent \textbf{Reformatted version of paper
    accepted for publication in:} \\ #4
   \vspace{15mm} \\ \noindent{\copyright} #5
   \clearpage
   \fancyhf{}
   \fancyhead[RE]{\sc{\nouppercase{#1}}}
   \fancyhead[LO]{\sc{\nouppercase{\rightmark}}}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% \makepapersubmitted command, same as above
% for submitted, not yet accepted paper
% ------------------------------------------
\newcommand{\makepaperacceptedat}[5]
{
   \fancyhead[RE]{}
   \fancyhead[LO]{}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
   \clearpage
   \thispagestyle{plain}
   \chapter[\space]{#2}
   \vspace{25mm}
   \textbf{Authors:}\\ #3
   \vspace{10mm} \\ \noindent \textbf{Reformatted version of paper
    accepted at:} \\ #4
   \vspace{15mm} \\ \noindent{\copyright} #5
   \clearpage
   \fancyhf{}
   \fancyhead[RE]{\sc{\nouppercase{#1}}}
   \fancyhead[LO]{\sc{\nouppercase{\rightmark}}}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
}


% ------------------------------------------
% \makepapersubmitted command, same as above
% for submitted, not yet accepted paper
% ------------------------------------------
\newcommand{\makepapersubmitted}[5]
{
   \fancyhead[RE]{}
   \fancyhead[LO]{}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
   \clearpage
   \thispagestyle{plain}
   \chapter[\space]{#2}
   \vspace{25mm}
   \textbf{Authors:}\\ #3
   \vspace{10mm} \\ \noindent \textbf{Reformatted version of paper
    submitted to:} \\ #4
   \vspace{15mm} \\ \noindent{\copyright} #5
   \clearpage
   \fancyhf{}
   \fancyhead[RE]{\sc{\nouppercase{#1}}}
   \fancyhead[LO]{\sc{\nouppercase{\rightmark}}}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% \makepapertobesubmitted command,
% for Manuscript, not yet submitted
% ------------------------------------------
\newcommand{\makepapertobesubmitted}[5]
{
   \fancyhead[RE]{}
   \fancyhead[LO]{}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
   \clearpage
   \thispagestyle{plain}
   \chapter[\space]{#2}
   \vspace{25mm}
   \textbf{Authors:}\\ #3
   \vspace{10mm} \\ \noindent \textbf{Manuscript to be submitted.} \\
   \clearpage
   \fancyhf{}
   \fancyhead[RE]{\sc{\nouppercase{#1}}}
   \fancyhead[LO]{\sc{\nouppercase{\rightmark}}}
   \fancyhead[LE]{\thepage} 
   \fancyhead[RO]{\thepage}
}

% ------------------------------------------
% \makebib command
% ------------------------------------------
\newcommand{\makebib}
{
    \def\bibname{}
    \clearpage
    \fancyhead[RE]{\sc References}
    \fancyhead[LO]{\sc References}
    \fancyhead[LE]{\thepage} 
    \fancyhead[RO]{\thepage}
    \thispagestyle{plain}
    \vspace*{2cm}
    \hrule
    \vspace*{8\p@}%
    \par
    {\Huge{\hfill\sc References}}\\[-6pt]
    \hrule
    \addcontentsline{toc}{chapter}{\sc References}
    \putbib
    \clearpage
}

% ------------------------------------------
% \makeappendices command
% ------------------------------------------
\newcommand{\makeappendices}
{
    \def\bibname{}
    \clearpage
    \fancyhead[RE]{\sc Appendices}
    \fancyhead[LO]{\sc Appendices}
    \fancyhead[LE]{\thepage} 
    \fancyhead[RO]{\thepage}
    \thispagestyle{plain}
    \vspace*{2cm}
    \hrule
    \vspace*{8\p@}%
    \par
    {\Huge{\hfill\sc Appendices}}\\[-6pt]
    \hrule
    \addcontentsline{toc}{chapter}{\sc Appendices}
    \clearpage
}

% ------------------------------------------
% \appendix command
% ------------------------------------------
\renewcommand{\appendix}%
{
  \clearpage%
  \setcounter{section}{0} 
  % ------------------------------------------
  % Adjust section headings
  % ------------------------------------------
  \renewcommand\section{\@startsection {section}{1}{-1pt}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\Large\bfseries}}

  \renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\large\bfseries}}
  \renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalsize\bfseries}}
  \renewcommand{\thesection}{\arabic{chapter}\Alph{section}}%
  \renewcommand{\thesubsection}{%
                \arabic{chapter}\Alph{section}.\arabic{subsection}}%
  \renewcommand{\theequation}{\arabic{chapter}\Alph{section}.\arabic{equation}}
}

% ------------------------------------------
% \paperappendix command
% ------------------------------------------
\newcommand{\paperappendix}%
{
  \setcounter{section}{0}%
  \renewcommand{\thesection}{\Alph{section}}%
  \renewcommand{\thesubsection}{\Alph{section}.\arabic{subsection}}%
  \renewcommand{\thesubsubsection}{%
                \Alph{section}.\arabic{subsection}.\arabic{\thesubsubsection}}%
  \renewcommand{\theequation}{\Alph{section}.\arabic{equation}}
}
