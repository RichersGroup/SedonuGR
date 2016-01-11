set grid
set term pdf
set output 'compare.pdf'
set yrange [-0.15:0.15]
#set log y
plot 'output.dat' every ::7 u 1:(($2-$4)/$4) title 'E relative error', \
     'output.dat' every ::7 u 1:(($3-$5)/$5) title 'N relative error', \
     'output.dat' every ::7 u 1:2 title 'E', \
     'output.dat' every ::7 u 1:3 title 'N'
     
