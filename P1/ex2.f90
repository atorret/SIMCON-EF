PROGRAM EX2
    IMPLICIT NONE
    REAL, PARAMETER :: pi = 3.141592653589
    REAL :: pi_aprox = 0.0000000, error = 1.0000000
    INTEGER :: K = 0
    INTEGER :: I
    REAL, DIMENSION(294337,2) :: results
    INCLUDE '../chdir_to_code.inc'

    DO WHILE (abs(error) > 0.001)
        pi_aprox = pi_aprox + (-1)**K * 4.0 / (2.0*K + 1)
        error = pi_aprox - pi
        K = K + 1
        results(K,1) = K
        results(K,2) = pi_aprox
    END DO
    OPEN(UNIT=10, FILE='ex2_results.dat', STATUS='REPLACE', ACTION='WRITE')
    DO I = 1, K
        WRITE(10, '(I8, 1X, ES16.8)') NINT(results(I,1)), results(I,2)
    END DO
    CLOSE(10)
    PRINT *, 'The approximate value of pi is for K =', K, 'is: ', pi_aprox
END PROGRAM EX2
