// Paper configuration
#let configuration(
 title: none,
 authors: (),
 abstract: [],
 cols: 2,
 enable_bibliography: false,
 enable_figure_outline: false,
 enable_table_outline: false,
 doc,
) = {
 // ----- Show rules -----
 // headings
 show heading: it => [
  #set align(left)
  #{
   if it.level > 1 {
    set text(10pt, weight: "semibold")
   } else {
    set text(12pt, weight: "bold")
   }
  }
  #block(smallcaps(it.body))
 ]
 // links and references
 show link: set text(blue)
 show ref: set text(blue)
 // math
 set math.equation(numbering: "(1)")
 show math.equation.where(block: false): box.with(
  fill: luma(250),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
 )
 show math.equation.where(block: true): block.with(
  width: 100%,
  fill: luma(250),
  inset: 8pt,
  radius: 3pt,
 )
 // code
 show raw.where(block: false): box.with(
  fill: luma(250),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 3pt),
  radius: 2pt,
 )
 show raw.where(block: true): block.with(
  width: 100%,
  fill: luma(250),
  inset: 8pt,
  radius: 3pt,
 )
 // figure
 show figure: set par(justify: false)
 show figure: set block(breakable: false)
 show figure.where(
  kind: image
 ): set image(width: 90%)

 // ------- Layout -------
 set page(
  paper: "a4",
  header: align(
   right + horizon,
   title,
  ),
  footer: [
   #set align(center)
   #context { counter(page).display("1", both: false) }
  ]
 )
 // diversions from default page layout
 set page(header: context {
  // don't show header on first page
  if here().page() > 1 {
   align(
    right + horizon,
    title,
   )
  } else {
   // show FAU logo on first page
   align(
    left,
    image("fau.png", width: 30%),
   )
  }
 })
 set text(lang: "en")
 set par(justify: true)
 set text(size: 10pt)

 // title
 set align(center)
 text(17pt, weight: "bold", title)

 // authors
 let count = authors.len()
 let ncols = calc.min(count, 3)
 grid(
  columns: (1fr,) * ncols,
  row-gutter: 24pt,
  ..authors.map(author => [
   #author.name \
   #link("mailto:" + author.email)
  ]),
 )

 // abstract
 set align(left)
 par(justify: false)[
  = Abstract
  #abstract
 ]
 
 // paper contents
 columns(cols, doc)

 if enable_bibliography or enable_figure_outline or enable_table_outline {
  pagebreak()
  // Bibliography
  if enable_bibliography {
   bibliography("bib.yml")
  }
  // List of figures
  if enable_figure_outline {
   outline(
    title: [List of figures],
    target: figure.where(kind: image),
   )
  }
  // List of tables
  if enable_table_outline {
   outline(
    title: [List of tables],
    target: figure.where(kind: table),
   )
  }
 }
}
