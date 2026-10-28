$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$output = Join-Path (Split-Path -Parent $here) 'output\pdf'
$job = 'DYN-Note_Experimentale_Grid_Peeling_v0_1'

New-Item -ItemType Directory -Force -Path $output | Out-Null
Push-Location $here
try {
    $latexArgs = @(
        '-interaction=nonstopmode',
        '-halt-on-error',
        '-file-line-error',
        '-synctex=0',
        "-jobname=$job",
        'main.tex'
    )

    & pdflatex @latexArgs
    if ($LASTEXITCODE -ne 0) { throw "Premier passage pdfLaTeX en échec ($LASTEXITCODE)." }

    & bibtex $job
    if ($LASTEXITCODE -ne 0) { throw "BibTeX en échec ($LASTEXITCODE)." }

    & pdflatex @latexArgs
    if ($LASTEXITCODE -ne 0) { throw "Deuxième passage pdfLaTeX en échec ($LASTEXITCODE)." }

    & pdflatex @latexArgs
    if ($LASTEXITCODE -ne 0) { throw "Troisième passage pdfLaTeX en échec ($LASTEXITCODE)." }

    Copy-Item -Force -LiteralPath "$job.pdf" -Destination (Join-Path $output "$job.pdf")
}
finally {
    Pop-Location
}
