BLOCK
    CHARACTER(LEN=512) :: exe
    INTEGER :: j, sep
    CALL GET_COMMAND_ARGUMENT(0, exe)
    sep = 0
    DO j = LEN_TRIM(exe), 1, -1
        IF (exe(j:j) == '/' .OR. exe(j:j) == '\') THEN
            sep = j
            EXIT
        END IF
    END DO
    IF (sep > 1) CALL CHDIR(exe(1:sep-1))
END BLOCK
