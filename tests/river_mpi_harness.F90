PROGRAM river_mpi_harness
   USE mpi, only: MPI_Init, MPI_Finalize, MPI_Comm_rank, MPI_Comm_size, &
      MPI_Comm_split, MPI_Comm_free, MPI_Barrier, MPI_Allreduce, &
      MPI_COMM_WORLD, MPI_COMM_NULL, MPI_UNDEFINED, MPI_INTEGER, MPI_SUM
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task, only: p_is_worker, p_comm_worker, p_iam_worker, &
      p_np_worker, p_np_io, p_err
   IMPLICIT NONE
   integer :: rank, nranks, local_failures, global_failures

   CALL MPI_Init(p_err)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, p_err)
   CALL MPI_Comm_size(MPI_COMM_WORLD, nranks, p_err)

   p_is_worker = .true.
   p_comm_worker = MPI_COMM_WORLD
   p_iam_worker = rank
   p_np_worker = nranks
   p_np_io = 0

   local_failures = 0
   CALL test_repeated_ordered_push(rank, nranks, local_failures)
   CALL test_zero_local_data_rank(rank, nranks, local_failures)
   CALL test_colm_role_split(rank, nranks, local_failures)
   CALL MPI_Allreduce(local_failures, global_failures, 1, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_WORLD, p_err)

   IF (rank == 0) THEN
      IF (global_failures == 0) THEN
         write(*,'(A,I0,A)') 'river MPI harness: PASS (', nranks, ' ranks)'
      ELSE
         write(*,'(A,I0)') 'river MPI harness: FAIL count=', global_failures
      ENDIF
   ENDIF
   CALL MPI_Finalize(p_err)
   IF (global_failures /= 0) ERROR STOP 1

CONTAINS

   SUBROUTINE test_colm_role_split(rank, nranks, failures)
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      integer :: color
      logical :: is_master, is_io

      ! Match the relevant CoLM topology: rank 0 is IO, the last rank is the
      ! master, and only intervening ranks participate in worker collectives.
      is_master = rank == nranks - 1
      is_io = rank == 0
      p_is_worker = .not. is_master .and. .not. is_io
      IF (p_is_worker) THEN
         color = 1
      ELSE
         color = MPI_UNDEFINED
      ENDIF
      CALL MPI_Comm_split(MPI_COMM_WORLD, color, rank, p_comm_worker, p_err)

      IF (p_is_worker) THEN
         CALL MPI_Comm_rank(p_comm_worker, p_iam_worker, p_err)
         CALL MPI_Comm_size(p_comm_worker, p_np_worker, p_err)
         CALL test_repeated_ordered_push(p_iam_worker, p_np_worker, failures)
         CALL test_zero_local_data_rank(p_iam_worker, p_np_worker, failures)
         CALL MPI_Comm_free(p_comm_worker, p_err)
      ELSE
         p_comm_worker = MPI_COMM_NULL
         p_iam_worker = -1
         p_np_worker = 0
         CALL test_nonworker_noop(failures)
      ENDIF

      ! In the two-rank case there are no workers.  This world barrier proves
      ! that IO/master ranks still leave the role-split path collectively.
      CALL MPI_Barrier(MPI_COMM_WORLD, p_err)
   END SUBROUTINE test_colm_role_split

   SUBROUTINE test_nonworker_noop(failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, build_worker_pushdata, worker_push_data
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      integer :: ids(1)
      real(r8) :: send_value(1), recv_value(1)

      ids = 1
      send_value = 7._r8
      recv_value = -77._r8
      CALL build_worker_pushdata(0, ids, 0, ids, push)
      CALL worker_push_data(push, send_value, recv_value, -1._r8)
      IF (recv_value(1) /= -77._r8) failures = failures + 1
   END SUBROUTINE test_nonworker_noop

   SUBROUTINE test_repeated_ordered_push(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push_next, push_previous
      integer :: ids_me(1), ids_next(1), ids_previous(1), iteration
      real(r8) :: send_value(1), recv_next(1), recv_previous(1)
      real(r8) :: expected_next, expected_previous

      ids_me(1) = rank + 1
      ids_next(1) = modulo(rank + 1, nranks) + 1
      ids_previous(1) = modulo(rank - 1 + nranks, nranks) + 1
      CALL build_worker_pushdata(1, ids_me, 1, ids_next, push_next)
      CALL build_worker_pushdata(1, ids_me, 1, ids_previous, push_previous)

      DO iteration = 1, 64
         send_value(1) = real(iteration * 1000 + ids_me(1), r8)
         recv_next = -1._r8
         recv_previous = -1._r8
         CALL worker_push_data(push_next, send_value, recv_next, -1._r8)
         CALL worker_push_data(push_previous, send_value, recv_previous, -1._r8)
         expected_next = real(iteration * 1000 + ids_next(1), r8)
         expected_previous = real(iteration * 1000 + ids_previous(1), r8)
         IF (recv_next(1) /= expected_next) failures = failures + 1
         IF (recv_previous(1) /= expected_previous) failures = failures + 1
      ENDDO
   END SUBROUTINE test_repeated_ordered_push

   SUBROUTINE test_zero_local_data_rank(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      integer :: num_me, request_id(1), iteration
      integer, allocatable :: ids_me(:)
      real(r8), allocatable :: send_value(:)
      real(r8) :: recv_value(1), expected

      IF (rank == 0) THEN
         num_me = 0
      ELSE
         num_me = 1
      ENDIF
      allocate(ids_me(num_me), send_value(num_me))
      IF (num_me == 1) ids_me(1) = rank
      request_id(1) = modulo(rank, nranks - 1) + 1
      CALL build_worker_pushdata(num_me, ids_me, 1, request_id, push)

      DO iteration = 1, 32
         IF (num_me == 1) send_value(1) = real(iteration * 1000 + ids_me(1), r8)
         recv_value = -1._r8
         CALL worker_push_data(push, send_value, recv_value, -1._r8)
         expected = real(iteration * 1000 + request_id(1), r8)
         IF (recv_value(1) /= expected) failures = failures + 1
      ENDDO
   END SUBROUTINE test_zero_local_data_rank
END PROGRAM river_mpi_harness
