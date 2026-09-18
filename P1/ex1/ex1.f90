PROGRAM EX1
    IMPLICIT NONE
    INTEGER :: A
    REAL    :: B, RESULT
    INCLUDE '../../chdir_to_code.inc'

    PRINT *, 'Enter an integer:'
    READ(*, *) A

    PRINT *, 'Enter a real number:'
    READ(*, *) B

    RESULT = REAL(A) + B

    PRINT '(A, F10.2)', 'The sum is: ', RESULT

    OPEN(UNIT=10, FILE='ex1_results.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE(10, '(A, F10.2)') 'The sum is: ', RESULT
    CLOSE(10)

END PROGRAM EX1