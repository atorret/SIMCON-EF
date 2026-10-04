PROGRAM EX1
    IMPLICIT NONE
    INTEGER :: A
    REAL    :: B, RESULT

    ! Sets the working directory to the source folder so programs can be run/compiled 
    ! directly from the root workspace in VS Code without manual navigation.
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