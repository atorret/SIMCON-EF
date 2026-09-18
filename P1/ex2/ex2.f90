PROGRAM EX2
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846D0
    DOUBLE PRECISION ::  pi_aprox = 0.0D0, error    = 1.0D0
    INTEGER ::K = 0, sign_K   = 1
    INCLUDE '../../chdir_to_code.inc'
    OPEN(UNIT=10, FILE='ex2_results.dat', STATUS='REPLACE', ACTION='WRITE')

    DO WHILE (ABS(error) > 1.0D-6)
        pi_aprox = pi_aprox + sign_K * 4.0D0 / (2.0D0 * DBLE(K) + 1.0D0)
        error = pi_aprox - pi
        
        K = K + 1
        WRITE(10, '(I8, 1X, ES16.8)') K, pi_aprox
        sign_K = -sign_K
    END DO

    CLOSE(10)

    PRINT *, 'The approximate value of pi for K =', K, 'is: ', pi_aprox
END PROGRAM EX2