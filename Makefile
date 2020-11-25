.PHONY: run

build: jl rinocchio

jl:
	g++ -O2 -std=c++11 src/main.cpp src/joy_libert.cpp -o out/main -pthread -lntl -lgmp -lm

rinocchio:
	g++ -O2 -std=c++11 src/rinocchio.cpp src/joy_libert.cpp -o out/rinocchio -pthread -lntl -lgmp -lm

run-setup:
	./out/rinocchio

test: rinocchio run-setup