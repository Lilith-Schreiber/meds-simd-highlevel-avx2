# MEDS

This is a separated README file that states the current states and error messages of MEDS digital signature implementation.

**Branch: main**
Using AVX512 intrinsics

**Branch: avx2**
Adapting AVX2 intrinsics

To run the code, please reference to README.md

To switch the branch, please use `git checkout avx2` to switch to AVX2 intrinsics version; `git checkout main` to switch to AVX512 intrinsics version

**Result of cycles**
Type `make RUN` under the directory of `/ref/build`, and the result would be printed as

```
keypair (normal): 1335145
keypair   (SIMD): 7836581

   sign (normal): 5111855
   sign   (SIMD): 6923394

 verify (normal): 4837717
 verify   (SIMD): 7717983

Time (min of 1 runs):
keygen: 0.001633   (7836581 cycles)
sign:   0.001442   (6923394 cycles)
verify: 0.001608   (7717983 cycles)
```

**Error messages**
The error happens in `verify`

```
ERROR: bistream - read esceeds buffer!
Signature verification failed!
!!! FAILED   (SIMD) !!!
```
