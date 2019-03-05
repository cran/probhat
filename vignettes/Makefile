default:
	make probhat.pdf
	make removeother

probhat.pdf: probhat.tex
	pdflatex probhat.tex

removeother:
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f temp*
	rm -f Rplots*
