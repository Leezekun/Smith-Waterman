#!/bin/bash

.PHONY: gpp
gpp:
	g++ SmithWaterman.cpp functions.cpp globals.cpp -o SmithWaterman -W -Wall -Wextra -fopenmp

.PHONY: gppO3
gppO3:
	g++ SmithWaterman.cpp functions.cpp globals.cpp -o SmithWaterman -W -Wall -Wextra -fopenmp -O3

.PHONY: icpc
icpc:
	icpc SmithWaterman.cpp functions.cpp globals.cpp -o SmithWaterman -W -Wall -Wextra -fopenmp

.PHONY: icpcO3
icpcO3:
	icpc SmithWaterman.cpp functions.cpp globals.cpp -o SmithWaterman -W -Wall -Wextra -fopenmp -O3
