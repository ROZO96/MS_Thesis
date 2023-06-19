#g++ main.cpp -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a
#./a.out

clang++ -Xpreprocessor -fopenmp -lomp -O3 -I /usr/local/opt/openblas/include /usr/local/opt/openblas/lib/libopenblas.a main3D.cpp -o lairds3D
./lairds3D

# clang++ -I /home/morton/local/include  -L /home/morton/local/lib -lopenblas -o lairds main3D.cpp

# cd output
# pythonw plot1D.py &

# cd output
# pythonw vor_plot2D.py

# cd ../../../exact_sod
# python sod_plot_rd.py &
