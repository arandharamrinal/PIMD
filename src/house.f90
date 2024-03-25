! Reduction of a real symmetric matrix, A, to tridiagonal
! form by the Householder method. This is followed by the
! evaluation of the eigenvalues and eigenvectors.
!
! A - The (N,N) matrix to be diagonalized; physical dimension (NP,NP).
! D - Array of dimension N containing eigenvalues on output; PD: (NP,NP).

SUBROUTINE HOUSEDIAG(A, N, NP, D)
    implicit none
    INTEGER, INTENT(IN) :: N, NP
    REAL(kind=8), DIMENSION(NP, NP), INTENT(INOUT) :: A
    REAL(kind=8), DIMENSION(NP), INTENT(OUT) :: D

    INTEGER :: I, J, K, L, M, ITER
    REAL(kind=8) :: SCALES, B,C, H, F, G, HH, DD, R, S, P
    REAL(kind=8), DIMENSION(NP) :: E

    IF (N > 1) THEN
        DO I = N, 2, -1
            L = I - 1
            H = 0.d0
            SCALES = 0.d0

            IF (L > 1) THEN
                DO K = 1, L
                    SCALES = SCALES + ABS(A(I, K))
                END DO

                IF (SCALES == 0.d0) THEN
                    E(I) = A(I, L)
                ELSE
                    DO K = 1, L
                        A(I, K) = A(I, K) / SCALES
                        H = H + A(I, K)**2
                    END DO

                    F = A(I, L)
                    G = -SIGN(SQRT(H), F)
                    E(I) = SCALES * G
                    H = H - F * G
                    A(I, L) = F - G
                    F = 0.d0

                    DO J = 1, L
                        A(J, I) = A(I, J) / H
                        G = 0.d0
                        DO K = 1, J
                            G = G + A(J, K) * A(I, K)
                        END DO

                        IF (L > J) THEN
                            DO K = J + 1, L
                                G = G + A(K, J) * A(I, K)
                            END DO
                        END IF

                        E(J) = G / H
                        F = F + E(J) * A(I, J)
                    END DO

                    HH = F / (H + H)

                    DO J = 1, L
                        F = A(I, J)
                        G = E(J) - HH * F
                        E(J) = G

                        DO K = 1, J
                            A(J, K) = A(J, K) - F * E(K) - G * A(I, K)
                        END DO
                    END DO
                END IF
            ELSE
                E(I) = A(I, L)
            END IF

            D(I) = H
        END DO
    END IF

    D(1) = 0.d0
    E(1) = 0.d0

    DO I = 1, N
        L = I - 1

        IF (D(I) /= 0.d0) THEN
            DO J = 1, L
                G = 0.d0
                DO K = 1, L
                    G = G + A(I, K) * A(K, J)
                END DO

                DO K = 1, L
                    A(K, J) = A(K, J) - G * A(K, I)
                END DO
            END DO
        END IF

        D(I) = A(I, I)
        A(I, I) = 1.d0

        IF (L >= 1) THEN
            DO J = 1, L
                A(I, J) = 0.d0
                A(J, I) = 0.d0
            END DO
        END IF
    END DO

    ! Now diagonalize the tridiagonal matrix produced above.
    IF (N > 1) THEN
        DO I = 2, N
            E(I-1) = E(I)
        END DO

        E(N) = 0.d0

        DO L = 1, N
            ITER = 0
1           DO M = L, N-1
                DD = ABS(D(M)) + ABS(D(M+1))
                IF (ABS(E(M)) + DD == DD) THEN
                    GO TO 2
                END IF
            END DO

            M = N
2           IF (M /= L) THEN
                IF (ITER == 30) THEN
                    PAUSE 'too many iterations'
                END IF

                ITER = ITER + 1
                G = (D(L+1) - D(L)) / (2.d0 * E(L))
                R = SQRT(G**2 + 1.d0)
                G = D(M) - D(L) + E(L) / (G + SIGN(R, G))
                S = 1.d0
                P = 0.d0

                DO I = M-1, L, -1
                    F = S * E(I)
                    B = C * E(I)

                    IF (ABS(F) >= ABS(G)) THEN
                        C = G / F
                        R = SQRT(C**2 + 1.d0)
                        E(I+1) = F * R
                        S = 1.d0 / R
                        C = C * S
                    ELSE
                        S = F / G
                        R = SQRT(S**2 + 1.d0)
                        E(I+1) = G * R
                        C = 1.d0 / R
                        S = S * C
                    END IF

                    G = D(I+1) - P
                    R = (D(I) - G) * S + 2.d0 * C * B
                    P = S * R
                    D(I+1) = G + P
                    G = C * R - B

                    DO K = 1, N
                        F = A(K, I+1)
                        A(K, I+1) = S * A(K, I) + C * F
                        A(K, I) = C * A(K, I) - S * F
                    END DO
                END DO

                D(L) = D(L) - P
                E(L) = G
                E(M) = 0.d0
                GO TO 1
            END IF
        END DO
    END IF

    ! Now sort the eigenvalues and eigenvectors in increasing order.
    DO I = 1, N-1
        K = I
        P = D(I)

        DO J = I+1, N
            IF (D(J) <= P) THEN
                K = J
                P = D(J)
            END IF
		END DO
        IF(K.NE.I)THEN
			D(K)=D(I)
			D(I)=P
			DO J=1,N
			  P=A(J,I)
			  A(J,I)=A(J,K)
			  A(J,K)=P
			END DO
        ENDIF
	ENDDO
END subroutine 

