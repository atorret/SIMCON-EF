PROGRAM EX5
    IMPLICIT NONE
    REAL :: m=0.200, k=2.0
    REAL :: t_max=10.0, dt, x_0=0.05, v_0=0.0, t_0=0.0
    REAL, ALLOCATABLE :: x_euler_001(:), v_euler_001(:), t_euler_001(:)
    REAL, ALLOCATABLE :: x_euler_01(:),  v_euler_01(:),  t_euler_01(:)
    REAL, ALLOCATABLE :: x_euler_02(:),  v_euler_02(:),  t_euler_02(:)
    REAL, ALLOCATABLE :: x_euler_predictor(:), v_euler_predictor(:), t_euler_predictor(:)
    REAL, ALLOCATABLE :: x_verlet(:), v_verlet(:), t_verlet(:)
    INTEGER :: n

    dt = 0.001
    n = INT(t_max/dt) + 1
    ALLOCATE(x_euler_001(n), v_euler_001(n), t_euler_001(n))
    CALL EULER(m, k, t_max, dt, x_0, v_0, t_0, x_euler_001, v_euler_001, t_euler_001)

    dt = 0.01
    n = INT(t_max/dt) + 1
    ALLOCATE(x_euler_01(n), v_euler_01(n), t_euler_01(n))
    CALL EULER(m, k, t_max, dt, x_0, v_0, t_0, x_euler_01, v_euler_01, t_euler_01)

    dt = 0.02
    n = INT(t_max/dt) + 1
    ALLOCATE(x_euler_02(n), v_euler_02(n), t_euler_02(n))
    ALLOCATE(x_euler_predictor(n), v_euler_predictor(n), t_euler_predictor(n))
    ALLOCATE(x_verlet(n), v_verlet(n), t_verlet(n))
    CALL EULER(m, k, t_max, dt, x_0, v_0, t_0, x_euler_02, v_euler_02, t_euler_02)
    CALL EULER_PREDICTOR(m, k, t_max, dt, x_0, v_0, t_0, x_euler_predictor, v_euler_predictor, t_euler_predictor)
    CALL VERLET(m, k, t_max, dt, x_0, v_0, t_0, x_verlet, v_verlet, t_verlet)

    CALL WRITE_TRAJ("euler_001.dat", t_euler_001, x_euler_001, v_euler_001, m, k)
    CALL WRITE_TRAJ("euler_01.dat",  t_euler_01,  x_euler_01,  v_euler_01,  m, k)
    CALL WRITE_TRAJ("euler_02.dat",  t_euler_02,  x_euler_02,  v_euler_02,  m, k)
    CALL WRITE_TRAJ("euler_predictor.dat", t_euler_predictor, x_euler_predictor, v_euler_predictor, m, k)
    CALL WRITE_TRAJ("verlet.dat", t_verlet, x_verlet, v_verlet, m, k)

    PRINT *, "Wrote euler_*.dat, euler_predictor.dat, verlet.dat"

CONTAINS

    SUBROUTINE EULER(m, k, t_max, dt, x_0, v_0, t_0, x, v, t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x, v, t
        REAL :: a
        INTEGER :: i

        x(1) = x_0
        v(1) = v_0
        t(1) = t_0

        DO i = 2, INT(t_max/dt) + 1
            a = -(k/m) * x(i-1)
            v(i) = v(i-1) + a * dt
            x(i) = x(i-1) + v(i-1) * dt + 0.5 * a * dt**2
            t(i) = t(i-1) + dt
        END DO
    END SUBROUTINE EULER

    SUBROUTINE VERLET(m, k, t_max, dt, x_0, v_0, t_0, x, v, t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x, v, t
        REAL :: a, a_new
        INTEGER :: i

        x(1) = x_0
        v(1) = v_0
        t(1) = t_0
        a = -(k/m) * x(1)

        DO i = 2, INT(t_max/dt) + 1
            x(i) = x(i-1) + v(i-1) * dt + 0.5 * a * dt**2
            a_new = -(k/m) * x(i)
            v(i) = v(i-1) + 0.5 * (a + a_new) * dt
            a = a_new
            t(i) = t(i-1) + dt
        END DO
        RETURN
    END SUBROUTINE VERLET

    SUBROUTINE EULER_PREDICTOR(m, k, t_max, dt, x_0, v_0, t_0, x, v, t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x, v, t
        REAL :: a, x_p, a_p
        INTEGER :: i

        x(1) = x_0
        v(1) = v_0
        t(1) = t_0

        DO i = 2, INT(t_max/dt) + 1
            a = -(k/m) * x(i-1)
            x_p = x(i-1) + v(i-1) * dt + 0.5 * a * dt**2
            a_p = -(k/m) * x_p
            a = 0.5 * (a_p + a)
            x(i) = x(i-1) + v(i-1) * dt + 0.5 * a * dt**2
            v(i) = v(i-1) + a * dt
            t(i) = t(i-1) + dt
        END DO
        RETURN
    END SUBROUTINE EULER_PREDICTOR

    SUBROUTINE WRITE_TRAJ(fname, t, x, v, m, k)
        CHARACTER(LEN=*), INTENT(IN) :: fname
        REAL, INTENT(IN) :: t(:), x(:), v(:), m, k
        REAL :: ek, ep, et
        INTEGER :: i, u
        u = 20
        OPEN(UNIT=u, FILE=fname, STATUS="REPLACE", ACTION="WRITE")
        WRITE(u, '(A)') "# t x v Ek Ep ET"
        DO i = 1, SIZE(t)
            ek = 0.5 * m * v(i)**2
            ep = 0.5 * k * x(i)**2
            et = ek + ep
            WRITE(u, '(6ES16.8)') t(i), x(i), v(i), ek, ep, et
        END DO
        CLOSE(u)
    END SUBROUTINE WRITE_TRAJ

END PROGRAM EX5