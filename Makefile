OUTDIR = out

build: joye_libert rinocchio

test: joye_libert
	./$(OUTDIR)/jltest

run: rinocchio
	./$(OUTDIR)/rinocchio

runinv: invertible
	./$(OUTDIR)/invertible

runint: interpolate
	./$(OUTDIR)/interpolate

joye_libert: | out
	g++ -O2 -std=c++11 -I./include test/joye_libert_test.cpp src/*.cpp -o $(OUTDIR)/jltest -pthread -lntl -lgmp -lm

rinocchio: | out
	g++ -O2 -std=c++11 -I./include test/rinocchio_test.cpp src/*.cpp -o $(OUTDIR)/rinocchio -pthread -lntl -lgmp -lm

invertible: | out
	g++ -O2 -std=c++11 -I./include test/invertible_test.cpp src/*.cpp -o $(OUTDIR)/invertible -pthread -lntl -lgmp -lm

interpolate: | out
	g++ -O2 -std=c++11 -I./include test/interpolate_test.cpp src/*.cpp -o $(OUTDIR)/interpolate -pthread -lntl -lgmp -lm

out:
	mkdir -p $(OUTDIR)
