.PHONY=clean

script: script.cpp 
	g++ -std=c++11 -Werror  script.cpp -o script

all_max: all_max.cpp
	g++ -std=c++11 -Werror  all_max.cpp -o all_max

test: test.cpp 
	g++ -std=c++11 -Werror  test.cpp -o test

clean:
	rm -f *.o script all_max test

