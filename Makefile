


.PHONY: all
all: animation.gif


.PHONY: clean
clean:
	rm -rf /tmp/anim
	rm -f animation.gif


main.out: main.cpp
	g++ $^ -o $@ -DNDEBUG -march=native -O3 -Wall -Wpedantic -std=c++17


animation.gif: main.out
	rm -rf /tmp/anim
	mkdir /tmp/anim
	./$< 0.25 1.5 0.9 2 100 100 500
	convert -delay 4 -loop 0 /tmp/anim/* $@
