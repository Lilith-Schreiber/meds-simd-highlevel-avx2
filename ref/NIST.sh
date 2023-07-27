#!/bin/bash

for p in `./params.py -l`; do
  v=Reference_Implementation

  mkdir -p MEDS/$v/$p/src
  mkdir -p MEDS/$v/$p/include

  cp src/*.c MEDS/$v/$p/src
  cp include/*.h MEDS/$v/$p/include

  rm MEDS/$v/$p/src/randombytes.c
  rm MEDS/$v/$p/src/KAT_test.c
  rm MEDS/$v/$p/include/randombytes.h

  cp src/NIST/PQCgenKAT_sign.c MEDS/$v/$p/src/PQCgenKAT_sign.c
  cp src/NIST/rng.c MEDS/$v/$p/src/randombytes.c
  cp include/NIST/rng.h MEDS/$v/$p/include/randombytes.h
  cp include/NIST/rng.h MEDS/$v/$p/include/rng.h

  ./params.py -p $p > MEDS/$v/$p/include/params.h

  ./params.py -a $p > MEDS/$v/$p/include/api.h

  cp NIST.cmk MEDS/$v/$p/CMakeLists.txt
done

ln -s $v MEDS/Optimized_Implementation

