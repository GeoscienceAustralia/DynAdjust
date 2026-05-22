# Static build optimizations
# This file contains common optimizations for static builds

function(optimize_static_target TARGET_NAME)
    # Enable LTO for this specific target if supported
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set_property(TARGET ${TARGET_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()

    # Platform-specific optimizations
    if(UNIX AND NOT APPLE)
        target_link_options(${TARGET_NAME} PRIVATE
            -static-libgcc
            -static-libstdc++
            -Wl,--as-needed
            -Wl,-O2
            -Wl,--strip-all
            # glibc 2.34+ merged libpthread into libc; libpthread.a is an
            # empty stub. libstdc++.a references ALL pthread functions via
            # weak aliases (__gthrw_), and glibc defines them as WEAK (W)
            # in libc.a. LLD will not pull a weak-defined archive member to
            # satisfy a weak-undefined reference, leaving every pthread
            # symbol at address 0 and crashing at runtime.
            # Fix: add a strong undefined reference (-u) for each symbol.
            # This forces LLD to pull in each object from libc.a, making
            # the weak-defined symbol available for the weak refs to bind to.
            -Wl,-u,pthread_once
            -Wl,-u,pthread_getspecific
            -Wl,-u,pthread_setspecific
            -Wl,-u,pthread_create
            -Wl,-u,pthread_join
            -Wl,-u,pthread_equal
            -Wl,-u,pthread_self
            -Wl,-u,pthread_detach
            -Wl,-u,pthread_cancel
            -Wl,-u,pthread_mutex_lock
            -Wl,-u,pthread_mutex_trylock
            -Wl,-u,pthread_mutex_unlock
            -Wl,-u,pthread_mutex_init
            -Wl,-u,pthread_mutex_destroy
            -Wl,-u,pthread_cond_init
            -Wl,-u,pthread_cond_broadcast
            -Wl,-u,pthread_cond_signal
            -Wl,-u,pthread_cond_wait
            -Wl,-u,pthread_cond_timedwait
            -Wl,-u,pthread_cond_destroy
            -Wl,-u,pthread_key_create
            -Wl,-u,pthread_key_delete
            -Wl,-u,pthread_mutexattr_init
            -Wl,-u,pthread_mutexattr_settype
            -Wl,-u,pthread_mutexattr_destroy
            -Wl,-u,pthread_attr_init
            -Wl,-u,pthread_attr_destroy
            -Wl,-u,pthread_attr_setdetachstate
            -Wl,-u,pthread_exit
            -Wl,-u,__pthread_key_create
        )
        target_link_libraries(${TARGET_NAME} PRIVATE pthread m dl)
        target_compile_options(${TARGET_NAME} PRIVATE
            -ffunction-sections
            -fdata-sections
        )
    elseif(APPLE)
        # macOS-specific optimizations
        target_link_options(${TARGET_NAME} PRIVATE
            -Wl,-dead_strip         # Remove unused code
            -Wl,-dead_strip_dylibs  # Remove unused dylibs
        )
    elseif(WIN32 AND MSVC)
        # Windows MSVC-specific optimizations
        target_compile_options(${TARGET_NAME} PRIVATE
            /Gy                     # Enable function-level linking
        )
        target_link_options(${TARGET_NAME} PRIVATE
            /OPT:REF                # Remove unreferenced code
            /OPT:ICF                # Identical COMDAT folding
        )
    endif()
endfunction()