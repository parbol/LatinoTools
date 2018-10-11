#include "SkimMeBaby.C"

//g++ -o runSkim mt2_bisect.cpp runSkim.C `root-config --cflags --libs`


int main(int argc, char **argv) {

    SkimMeBaby *a = new SkimMeBaby(argv[1], argv[2]);
    return 1;

}
