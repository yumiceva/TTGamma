all: ttgamma

ttgamma: ttgamma3Dict.o ttgamma3.o runttgamma.C
	g++ -o ttgamma `root-config --libs` -I`root-config --incdir` ttgamma3Dict.o ttgamma3.o runttgamma.C

EventTree.o: EventTree.cpp EventTree.h
	g++ -c -I`root-config --incdir` EventTree.cpp

Selector.o: Selector.cpp Selector.h
	g++ -c -I`root-config --incdir` Selector.cpp

ttgamma3.o: ttgamma3.C ttgamma3.h
	g++ -c -I`root-config --incdir` ttgamma3.C

ttgamma3Dict.o: ttgamma3Dict.C ttgamma3Dict.h
	g++ -c -I`root-config --incdir` ttgamma3Dict.C

ttgamma3Dict.C: ttgamma3.h
	rootcint -f ttgamma3Dict.C -c ttgamma3.h Linkdef.h

clean:
	rm ttgamma3.o ttgamma3Dict.o ttgamma3Dict.h ttgamma3Dict.C
