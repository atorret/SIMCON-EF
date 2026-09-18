PROGRAM EX5
    IMPLICIT NONE
    REAL :: m=0.200,k=2.0
    REAL :: t_max=10.0, dt, x_0=0.05, v_0=0.0, t_0=0.0
    REAL, ALLOCATABLE :: x_euler_001(:),v_euler_001(:), t_euler_001(:)
    REAL, ALLOCATABLE :: x_euler_01(:),v_euler_01(:), t_euler_01(:)
    REAL, ALLOCATABLE :: x_euler_02(:),v_euler_02(:), t_euler_02(:)

    REAL, ALLOCATABLE :: x_euler_predictor(:),v_euler_predictor(:), t_euler_predictor(:)
    REAL, ALLOCATABLE :: x_verlet(:),v_verlet(:), t_verlet(:)
    INTEGER :: i

    dt=0.001
    CALL EULER(m,k,t_max,dt,x_0,v_0,t_0,x_euler_001,v_euler_001,t_euler_001)
    OPEN(10,file='euler_001.dat')
    WRITE(10,*) 't,x,v'
    DO i = 1, INT(t_max/dt)+1
        WRITE(10,*) t_euler_001(i),x_euler_001(i),v_euler_001(i)
    END DO
    CLOSE(10)

    dt=0.01
    CALL EULER(m,k,t_max,dt,x_0,v_0,t_0,x_euler_01,v_euler_01,t_euler_01)
    OPEN(11,file='euler_01.dat')
    WRITE(11,*) 't,x,v'
    DO i = 1, INT(t_max/dt)+1
        WRITE(11,*) t_euler_01(i),x_euler_01(i),v_euler_01(i)
    END DO
    CLOSE(11)

    dt=0.02
    CALL EULER(m,k,t_max,dt,x_0,v_0,t_0,x_euler_02,v_euler_02,t_euler_02)
    OPEN(12,file='euler_02.dat')
    WRITE(12,*) 't,x,v'
    DO i = 1, INT(t_max/dt)+1
        WRITE(12,*) t_euler_02(i),x_euler_02(i),v_euler_02(i)
    END DO
    CLOSE(12)

    CALL EULER_PREDICTOR(m,k,t_max,dt,x_0,v_0,t_0,x_euler_predictor,v_euler_predictor,t_euler_predictor)
    OPEN(13,file='euler_predictor.dat')
    WRITE(13,*) 't,x,v'
    DO i = 1, INT(t_max/dt)+1
        WRITE(13,*) t_euler_predictor(i),x_euler_predictor(i),v_euler_predictor(i)
    END DO
    CLOSE(13)

    CALL VERLET(m,k,t_max,dt,x_0,v_0,t_0,x_verlet,v_verlet,t_verlet)
    OPEN(14,file='verlet.dat')
    WRITE(14,*) 't,x,v'
    DO i = 1, INT(t_max/dt)+1
        WRITE(14,*) t_verlet(i),x_verlet(i),v_verlet(i)
    END DO
    CLOSE(14)

CONTAINS
    SUBROUTINE EULER(m,k,t_max,dt,x_0,v_0,t_0,x,v,t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m,k,t_max,dt,x_0,v_0,t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x,v,t
        REAL :: a
        INTEGER :: i

        x(1)=x_0
        v(1)=v_0
        t(1)=t_0

        DO i = 2, INT(t_max/dt)+1
            a=-(k/m)*x(i-1)
            v(i)=v(i-1)+a*dt
            x(i)=x(i-1)+v(i-1)*dt+0.5*a*dt**2
            t(i)=t(i-1)+dt
        END DO
    END SUBROUTINE EULER

    SUBROUTINE VERLET(m, k, t_max, dt, x_0, v_0, t_0, x, v, t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m, k, t_max, dt, x_0, v_0, t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x, v, t
        REAL :: a
        INTEGER :: i, n
    
        n = NINT((t_max - t_0) / dt) + 1
    
        t(1) = t_0
        x(1) = x_0
        v(1) = v_0
        a    = -(k / m) * x(1)
    
        t(2) = t_0 + dt
        x(2) = x(1) + v(1) * dt + 0.5 * a * (dt**2)
        v(2) = v(1) + a * dt
    
        DO i = 3, n
            a    = -(k / m) * x(i-1)
            x(i) = 2.0 * x(i-1) - x(i-2) + a * (dt**2)
            t(i) = t(i-1) + dt
            v(i-1) = (x(i) - x(i-2)) / (2.0 * dt)
        END DO
    
        v(n) = (x(n) - x(n - 1)) / dt
    
        RETURN
    END SUBROUTINE VERLET

    SUBROUTINE EULER_PREDICTOR(m,k,t_max,dt,x_0,v_0,t_0,x,v,t)
        IMPLICIT NONE
        REAL, INTENT(IN) :: m,k,t_max,dt,x_0,v_0,t_0
        REAL, DIMENSION(:), INTENT(OUT) :: x,v,t
        REAL :: a, x_p, a_p
        INTEGER :: i

        x(1)=x_0
        v(1)=v_0
        t(1)=t_0

        DO i = 2, NINT(t_max/dt)+1
            a=-(k/m)*x(i-1)
            x_p=x(i-1)+v(i-1)*dt+0.5*a*dt**2
            a_p=-(k/m)*x_p
            a=0.5*(a_p+a)
            x(i)=x(i-1)+v(i-1)*dt+0.5*a*dt**2
            v(i)=v(i-1)+a*dt
            t(i)=t(i-1)+dt
        END DO

        RETURN
    END SUBROUTINE EULER_PREDICTOR

END PROGRAM EX5