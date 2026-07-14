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
   CALL test_single_batch_matches_scalar(rank, nranks, local_failures)
   CALL test_single_batch_zero_local_data(rank, nranks, local_failures)
   CALL test_batch_matches_scalar(rank, nranks, local_failures)
   CALL test_batch_zero_local_data(rank, nranks, local_failures)
   CALL test_multi_average_ignores_nonpositive_area(rank, nranks, local_failures)
   CALL test_multi_zero_source_slots(rank, local_failures)
   CALL test_remap_average_ignores_nonpositive_area(local_failures)
   CALL test_true_empty_worker_mapping(rank, local_failures)
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

   SUBROUTINE test_single_batch_matches_scalar(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      type(worker_push_real8_field_type) :: fields(3)
      integer :: ids_me(1), ids_req(1), iteration
      real(r8), target :: send_1(1), send_2(1), send_3(1)
      real(r8), target :: recv_1(1), recv_2(1), recv_3(1)
      real(r8) :: scalar_1(1), scalar_2(1), scalar_3(1)

      ids_me(1) = rank + 1
      ids_req(1) = modulo(rank + 1, nranks) + 1
      CALL build_worker_pushdata(1, ids_me, 1, ids_req, push)

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -71._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -72._r8
      fields(3)%send => send_3
      fields(3)%recv => recv_3
      fields(3)%fillvalue = -73._r8

      DO iteration = 1, 32
         send_1(1) = real(iteration * 1000 + ids_me(1), r8)
         send_2(1) = real(iteration * 2000 + 2 * ids_me(1), r8)
         send_3(1) = real(iteration * 3000 + 3 * ids_me(1), r8)

         CALL worker_push_data(push, send_1, scalar_1, fields(1)%fillvalue)
         CALL worker_push_data(push, send_2, scalar_2, fields(2)%fillvalue)
         CALL worker_push_data(push, send_3, scalar_3, fields(3)%fillvalue)

         recv_1 = -301._r8
         recv_2 = -302._r8
         recv_3 = -303._r8
         CALL worker_push_data(push, fields)

         IF (recv_1(1) /= scalar_1(1)) failures = failures + 1
         IF (recv_2(1) /= scalar_2(1)) failures = failures + 1
         IF (recv_3(1) /= scalar_3(1)) failures = failures + 1
      ENDDO
   END SUBROUTINE test_single_batch_matches_scalar

   SUBROUTINE test_single_batch_zero_local_data(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      type(worker_push_real8_field_type) :: fields(2)
      integer :: num_me, request_id(1), iteration
      integer, allocatable :: ids_me(:)
      real(r8), allocatable, target :: send_1(:), send_2(:)
      real(r8), target :: recv_1(1), recv_2(1)
      real(r8) :: scalar_1(1), scalar_2(1)

      IF (rank == 0) THEN
         num_me = 0
      ELSE
         num_me = 1
      ENDIF
      allocate(ids_me(num_me), send_1(num_me), send_2(num_me))
      IF (num_me == 1) ids_me(1) = rank
      request_id(1) = modulo(rank, nranks - 1) + 1
      CALL build_worker_pushdata(num_me, ids_me, 1, request_id, push)

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -81._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -82._r8

      DO iteration = 1, 16
         IF (num_me == 1) THEN
            send_1(1) = real(iteration * 100 + ids_me(1), r8)
            send_2(1) = real(iteration * 200 + ids_me(1), r8)
         ENDIF

         CALL worker_push_data(push, send_1, scalar_1, fields(1)%fillvalue)
         CALL worker_push_data(push, send_2, scalar_2, fields(2)%fillvalue)

         recv_1 = -401._r8
         recv_2 = -402._r8
         CALL worker_push_data(push, fields)

         IF (recv_1(1) /= scalar_1(1)) failures = failures + 1
         IF (recv_2(1) /= scalar_2(1)) failures = failures + 1
      ENDDO
   END SUBROUTINE test_single_batch_zero_local_data

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

   SUBROUTINE test_batch_matches_scalar(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      type(worker_push_real8_field_type) :: fields(3)
      integer :: ids_me(1), ids_req(3,1), iteration
      real(r8) :: area_req(3,1)
      real(r8), target :: send_1(1), send_2(1), send_3(1)
      real(r8), target :: recv_1(1), recv_2(1), recv_3(1)
      real(r8) :: scalar_1(1), scalar_2(1), scalar_3(1)

      ids_me(1) = rank + 1
      ids_req(1,1) = ids_me(1)
      ids_req(2,1) = modulo(rank - 1 + nranks, nranks) + 1
      ids_req(3,1) = ids_req(2,1)
      area_req(:,1) = [1._r8, 2._r8, 3._r8]
      CALL build_worker_pushdata(1, ids_me, 1, ids_req, area_req, push)

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -11._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -22._r8
      fields(3)%send => send_3
      fields(3)%recv => recv_3
      fields(3)%fillvalue = -33._r8

      DO iteration = 1, 32
         send_1(1) = real(iteration * 1000 + ids_me(1), r8)
         send_2(1) = real(iteration * 2000 + 2 * ids_me(1), r8)
         send_3(1) = real(iteration * 3000 + 3 * ids_me(1), r8)

         CALL worker_push_data(push, send_1, scalar_1, fields(1)%fillvalue, mode='sum')
         CALL worker_push_data(push, send_2, scalar_2, fields(2)%fillvalue, mode='sum')
         CALL worker_push_data(push, send_3, scalar_3, fields(3)%fillvalue, mode='sum')

         recv_1 = -101._r8
         recv_2 = -102._r8
         recv_3 = -103._r8
         CALL worker_push_data(push, fields, mode='sum')

         IF (recv_1(1) /= scalar_1(1)) failures = failures + 1
         IF (recv_2(1) /= scalar_2(1)) failures = failures + 1
         IF (recv_3(1) /= scalar_3(1)) failures = failures + 1
      ENDDO
   END SUBROUTINE test_batch_matches_scalar

   SUBROUTINE test_batch_zero_local_data(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      type(worker_push_real8_field_type) :: fields(3)
      integer :: num_me, request_id(1,1), iteration
      integer, allocatable :: ids_me(:)
      real(r8) :: area_req(1,1)
      real(r8), allocatable, target :: send_1(:), send_2(:), send_3(:)
      real(r8), target :: recv_1(1), recv_2(1), recv_3(1)
      real(r8) :: scalar_1(1), scalar_2(1), scalar_3(1)

      IF (rank == 0) THEN
         num_me = 0
      ELSE
         num_me = 1
      ENDIF
      allocate(ids_me(num_me), send_1(num_me), send_2(num_me), send_3(num_me))
      IF (num_me == 1) ids_me(1) = rank
      request_id(1,1) = modulo(rank, nranks - 1) + 1
      area_req(1,1) = 1._r8
      CALL build_worker_pushdata(num_me, ids_me, 1, request_id, area_req, push)

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -41._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -42._r8
      fields(3)%send => send_3
      fields(3)%recv => recv_3
      fields(3)%fillvalue = -43._r8

      DO iteration = 1, 16
         IF (num_me == 1) THEN
            send_1(1) = real(iteration * 100 + ids_me(1), r8)
            send_2(1) = real(iteration * 200 + ids_me(1), r8)
            send_3(1) = real(iteration * 300 + ids_me(1), r8)
         ENDIF

         CALL worker_push_data(push, send_1, scalar_1, fields(1)%fillvalue, mode='sum')
         CALL worker_push_data(push, send_2, scalar_2, fields(2)%fillvalue, mode='sum')
         CALL worker_push_data(push, send_3, scalar_3, fields(3)%fillvalue, mode='sum')

         recv_1 = -201._r8
         recv_2 = -202._r8
         recv_3 = -203._r8
         CALL worker_push_data(push, fields, mode='sum')

         IF (recv_1(1) /= scalar_1(1)) failures = failures + 1
         IF (recv_2(1) /= scalar_2(1)) failures = failures + 1
         IF (recv_3(1) /= scalar_3(1)) failures = failures + 1
      ENDDO
   END SUBROUTINE test_batch_zero_local_data

   SUBROUTINE test_multi_average_ignores_nonpositive_area(rank, nranks, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank, nranks
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push
      type(worker_push_real8_field_type) :: fields(2)
      integer :: ids_me(1), ids_req(4,2), previous_id
      real(r8) :: area_req(4,2), scalar_1(2), scalar_2(2)
      real(r8), target :: send_1(1), send_2(1), recv_1(2), recv_2(2)
      real(r8) :: previous_1, previous_2, expected_1, expected_2

      ids_me(1) = rank + 1
      previous_id = modulo(rank - 1 + nranks, nranks) + 1
      ids_req(:,1) = [ids_me(1), previous_id, previous_id, ids_me(1)]
      ids_req(:,2) = ids_req(:,1)
      area_req(:,1) = [0._r8, -1._r8, 0._r8, -2._r8]
      area_req(:,2) = [2._r8, 3._r8, 0._r8, -1._r8]
      CALL build_worker_pushdata(1, ids_me, 2, ids_req, area_req, push)

      send_1(1) = real(100 + ids_me(1), r8)
      send_2(1) = real(1000 + 10 * ids_me(1), r8)
      previous_1 = real(100 + previous_id, r8)
      previous_2 = real(1000 + 10 * previous_id, r8)
      expected_1 = (2._r8 * send_1(1) + 3._r8 * previous_1) / 5._r8
      expected_2 = (2._r8 * send_2(1) + 3._r8 * previous_2) / 5._r8

      CALL worker_push_data(push, send_1, scalar_1, -71._r8, mode='average')
      CALL worker_push_data(push, send_2, scalar_2, -72._r8, mode='average')

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -71._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -72._r8
      CALL worker_push_data(push, fields, mode='average')

      IF (scalar_1(1) /= fields(1)%fillvalue) failures = failures + 1
      IF (scalar_2(1) /= fields(2)%fillvalue) failures = failures + 1
      IF (recv_1(1) /= fields(1)%fillvalue) failures = failures + 1
      IF (recv_2(1) /= fields(2)%fillvalue) failures = failures + 1
      IF (abs(scalar_1(2)-expected_1) > 1.e-12_r8) failures = failures + 1
      IF (abs(scalar_2(2)-expected_2) > 1.e-12_r8) failures = failures + 1
      IF (abs(recv_1(2)-expected_1) > 1.e-12_r8) failures = failures + 1
      IF (abs(recv_2(2)-expected_2) > 1.e-12_r8) failures = failures + 1
   END SUBROUTINE test_multi_average_ignores_nonpositive_area

   SUBROUTINE test_remap_average_ignores_nonpositive_area(failures)
      USE MOD_WorkerPushData, only: worker_remapdata_type, worker_remap_data_grid2pset
      integer, intent(inout) :: failures
      type(worker_remapdata_type) :: remap
      real(r8) :: input(3), output(2), expected

      remap%npset = 2
      allocate(remap%npart(2), remap%part_to(2), remap%areapart(2))
      remap%npart = [2, 3]
      allocate(remap%part_to(1)%val(2), remap%areapart(1)%val(2))
      allocate(remap%part_to(2)%val(3), remap%areapart(2)%val(3))
      remap%part_to(1)%val = [1, 2]
      remap%areapart(1)%val = [0._r8, 0._r8]
      remap%part_to(2)%val = [1, 2, 3]
      remap%areapart(2)%val = [2._r8, 3._r8, -4._r8]

      input = [10._r8, 20._r8, 30._r8]
      output = -999._r8
      CALL worker_remap_data_grid2pset(remap, input, output, -77._r8, 'average')

      expected = (2._r8 * input(1) + 3._r8 * input(2)) / 5._r8
      IF (output(1) /= -77._r8) failures = failures + 1
      IF (abs(output(2)-expected) > 1.e-12_r8) failures = failures + 1
   END SUBROUTINE test_remap_average_ignores_nonpositive_area

   SUBROUTINE test_multi_zero_source_slots(rank, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank
      integer, intent(inout) :: failures
      integer, parameter :: nrecv = 3
      type(worker_pushdata_type) :: push
      type(worker_push_real8_field_type) :: fields(2)
      integer :: ids_me(1), ids_req(0,nrecv)
      real(r8) :: area_req(0,nrecv)
      real(r8), target :: send_1(1), send_2(1)
      real(r8), target :: recv_1(nrecv), recv_2(nrecv)
      real(r8) :: scalar_1(nrecv), scalar_2(nrecv)

      ! A valid no-upstream mapping has output columns but no source slots.
      ! Both APIs must overwrite every output with its field-specific fill.
      ids_me(1) = rank + 1
      send_1(1) = real(100 + rank, r8)
      send_2(1) = real(200 + rank, r8)
      CALL build_worker_pushdata(1, ids_me, nrecv, ids_req, area_req, push)

      scalar_1 = -901._r8
      scalar_2 = -902._r8
      CALL worker_push_data(push, send_1, scalar_1, -91._r8, mode='sum')
      CALL worker_push_data(push, send_2, scalar_2, -92._r8, mode='sum')

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -91._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -92._r8
      recv_1 = -911._r8
      recv_2 = -912._r8
      CALL worker_push_data(push, fields, mode='sum')

      IF (any(scalar_1 /= fields(1)%fillvalue)) failures = failures + 1
      IF (any(scalar_2 /= fields(2)%fillvalue)) failures = failures + 1
      IF (any(recv_1 /= fields(1)%fillvalue)) failures = failures + 1
      IF (any(recv_2 /= fields(2)%fillvalue)) failures = failures + 1
   END SUBROUTINE test_multi_zero_source_slots

   SUBROUTINE test_true_empty_worker_mapping(rank, failures)
      USE MOD_WorkerPushData, only: worker_pushdata_type, worker_push_real8_field_type, &
         build_worker_pushdata, worker_push_data
      integer, intent(in) :: rank
      integer, intent(inout) :: failures
      type(worker_pushdata_type) :: push_single, push_multi
      type(worker_push_real8_field_type) :: fields(3)
      integer :: num_me, num_req
      integer, allocatable :: ids_me(:), ids_req_single(:), ids_req_multi(:,:)
      real(r8), allocatable :: area_req(:,:)
      real(r8), allocatable, target :: send_1(:), send_2(:), send_3(:)
      real(r8), allocatable, target :: recv_1(:), recv_2(:), recv_3(:)
      real(r8), allocatable :: scalar_1(:), scalar_2(:), scalar_3(:), recv_single(:)

      ! Reproduce a Flow-shaped empty worker: rank 0 owns no data and requests
      ! no output, while the other ranks still build and use the same mapping
      ! collectively.  Both single and multi builders must accept true 0:0.
      IF (rank == 0) THEN
         num_me = 0
         num_req = 0
      ELSE
         num_me = 1
         num_req = 1
      ENDIF

      allocate(ids_me(num_me), ids_req_single(num_req), ids_req_multi(1,num_req))
      allocate(area_req(1,num_req))
      allocate(send_1(num_me), send_2(num_me), send_3(num_me))
      allocate(recv_1(num_req), recv_2(num_req), recv_3(num_req))
      allocate(scalar_1(num_req), scalar_2(num_req), scalar_3(num_req), recv_single(num_req))

      IF (num_me > 0) THEN
         ids_me(1) = rank
         send_1(1) = real(100 + rank, r8)
         send_2(1) = real(200 + rank, r8)
         send_3(1) = real(300 + rank, r8)
      ENDIF
      IF (num_req > 0) THEN
         ids_req_single(1) = rank
         ids_req_multi(1,1) = rank
         area_req(1,1) = 1._r8
      ENDIF

      CALL build_worker_pushdata(num_me, ids_me, num_req, ids_req_single, push_single)
      CALL worker_push_data(push_single, send_1, recv_single, -51._r8)
      IF (num_req > 0) THEN
         IF (recv_single(1) /= send_1(1)) failures = failures + 1
      ENDIF

      fields(1)%send => send_1
      fields(1)%recv => recv_1
      fields(1)%fillvalue = -61._r8
      fields(2)%send => send_2
      fields(2)%recv => recv_2
      fields(2)%fillvalue = -62._r8
      fields(3)%send => send_3
      fields(3)%recv => recv_3
      fields(3)%fillvalue = -63._r8
      CALL worker_push_data(push_single, send_1, scalar_1, fields(1)%fillvalue)
      CALL worker_push_data(push_single, send_2, scalar_2, fields(2)%fillvalue)
      CALL worker_push_data(push_single, send_3, scalar_3, fields(3)%fillvalue)
      CALL worker_push_data(push_single, fields)
      IF (any(recv_1 /= scalar_1)) failures = failures + 1
      IF (any(recv_2 /= scalar_2)) failures = failures + 1
      IF (any(recv_3 /= scalar_3)) failures = failures + 1

      CALL build_worker_pushdata(num_me, ids_me, num_req, ids_req_multi, area_req, push_multi)

      CALL worker_push_data(push_multi, send_1, scalar_1, fields(1)%fillvalue, mode='sum')
      CALL worker_push_data(push_multi, send_2, scalar_2, fields(2)%fillvalue, mode='sum')
      CALL worker_push_data(push_multi, send_3, scalar_3, fields(3)%fillvalue, mode='sum')
      CALL worker_push_data(push_multi, fields, mode='sum')

      IF (any(recv_1 /= scalar_1)) failures = failures + 1
      IF (any(recv_2 /= scalar_2)) failures = failures + 1
      IF (any(recv_3 /= scalar_3)) failures = failures + 1
   END SUBROUTINE test_true_empty_worker_mapping

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
