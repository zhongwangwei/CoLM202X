PROGRAM worker_push_mapping_benchmark
   USE mpi, only: MPI_Init, MPI_Finalize, MPI_Comm_rank, MPI_Comm_size, MPI_Barrier, &
      MPI_Allgather, MPI_Allgatherv, MPI_Alltoall, MPI_Alltoallv, MPI_Allreduce, &
      MPI_Wtime, MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MAX, MPI_SUM
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task, only: p_is_worker, p_comm_worker, p_iam_worker, p_np_worker, p_np_io, p_err
   USE MOD_Utils, only: insert_into_sorted_list1, find_in_sorted_list1, quicksort
   IMPLICIT NONE

   integer, parameter :: timing_trials = 3
   ! The production enum is private, while the mapping components are public.
   ! These values are used only to exercise the production push path with a
   ! test-built, semantically equivalent descriptor.
   integer, parameter :: mapping_single = 1, mapping_multi = 2
   integer :: rank, nranks, benchmark_failures

   CALL MPI_Init(p_err)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, p_err)
   CALL MPI_Comm_size(MPI_COMM_WORLD, nranks, p_err)

   p_is_worker = .true.
   p_comm_worker = MPI_COMM_WORLD
   p_iam_worker = rank
   p_np_worker = nranks
   p_np_io = 0
   benchmark_failures = 0

   CALL warmup_mapping()

   IF (rank == 0) THEN
      write(*,'(A)') 'worker_push_mapping_benchmark_v3'
      write(*,'(A,I0)') 'timing_trials=', timing_trials
      write(*,'(A)') 'implementation,kind,ranks,nlocal,nrequest,build_repeats,push_repeats,' // &
         'build_seconds_median,push_seconds_median,directory_temp_bytes_max,failures'
   ENDIF

   CALL run_suite(.false.)
   CALL run_suite(.true.)

   CALL MPI_Finalize(p_err)
   IF (benchmark_failures /= 0) ERROR STOP 1

CONTAINS

   SUBROUTINE run_suite(use_directory)
      logical, intent(in) :: use_directory

      CALL run_single_case('single-zero', 0, 0, 64, 256, .false., use_directory)
      CALL run_single_case('single-empty-rank', 64, 128, 24, 256, .true., use_directory)
      CALL run_single_case('single-duplicate', 64, 256, 24, 256, .false., use_directory)
      CALL run_single_case('single-small', 64, 128, 24, 256, .false., use_directory)
      CALL run_single_case('single-medium', 512, 1024, 8, 128, .false., use_directory)
      CALL run_single_case('single-large', 4096, 8192, 3, 32, .false., use_directory)

      CALL run_multi_case('multi-zero', 0, 0, 64, 256, .false., use_directory)
      CALL run_multi_case('multi-empty-rank', 64, 64, 24, 256, .true., use_directory)
      CALL run_multi_case('multi-duplicate', 64, 128, 24, 256, .false., use_directory)
      CALL run_multi_case('multi-small', 64, 64, 24, 256, .false., use_directory)
      CALL run_multi_case('multi-medium', 512, 512, 8, 128, .false., use_directory)
      CALL run_multi_case('multi-large', 4096, 4096, 3, 32, .false., use_directory)
   END SUBROUTINE run_suite

   SUBROUTINE warmup_mapping()
      USE MOD_WorkerPushData, only: worker_pushdata_type, build_worker_pushdata, worker_push_data
      type(worker_pushdata_type) :: ring_push, directory_push
      integer :: ids_me(1), ids_req(1)
      integer :: temp_bytes
      real(r8) :: send(1), recv(1)

      ids_me(1) = rank + 1
      ids_req(1) = modulo(rank + 1, nranks) + 1
      send(1) = real(ids_me(1), r8)
      CALL build_worker_pushdata(1, ids_me, 1, ids_req, ring_push)
      CALL worker_push_data(ring_push, send, recv, -1._r8)
      CALL build_directory_single(1, ids_me, 1, ids_req, directory_push, temp_bytes)
      CALL worker_push_data(directory_push, send, recv, -1._r8)
      CALL MPI_Barrier(MPI_COMM_WORLD, p_err)
   END SUBROUTINE warmup_mapping

   SUBROUTINE run_single_case(label, base_nlocal, nrequest, build_repeats, push_repeats, empty_rank, use_directory)
      USE MOD_WorkerPushData, only: worker_pushdata_type, build_worker_pushdata, worker_push_data
      character(len=*), intent(in) :: label
      integer, intent(in) :: base_nlocal, nrequest, build_repeats, push_repeats
      logical, intent(in) :: empty_rank, use_directory

      type(worker_pushdata_type), allocatable :: push
      integer, allocatable :: ids_me(:), ids_req(:)
      real(r8), allocatable :: send(:), recv(:)
      integer :: nlocal, iteration, trial, failures, failures_global
      integer :: build_repeats_used, push_repeats_used
      integer :: directory_temp_bytes, directory_temp_bytes_global
      real(r8) :: t0, elapsed, elapsed_global
      real(r8) :: build_samples(timing_trials), push_samples(timing_trials)
      real(r8) :: build_time, push_time

      CALL build_case_ids(base_nlocal, nrequest, empty_rank, nlocal, ids_me, ids_req)
      allocate(send(nlocal), recv(nrequest))
      IF (nlocal > 0) send = real(ids_me, r8)
      CALL benchmark_repeat_counts(build_repeats, push_repeats, build_repeats_used, push_repeats_used)

      DO trial = 1, timing_trials
         CALL MPI_Barrier(MPI_COMM_WORLD, p_err)
         t0 = MPI_Wtime()
         DO iteration = 1, build_repeats_used
            allocate(push)
            IF (use_directory) THEN
               CALL build_directory_single(nlocal, ids_me, nrequest, ids_req, push, directory_temp_bytes)
            ELSE
               CALL build_worker_pushdata(nlocal, ids_me, nrequest, ids_req, push)
            ENDIF
            deallocate(push)
         ENDDO
         elapsed = (MPI_Wtime() - t0) / real(build_repeats_used, r8)
         CALL MPI_Allreduce(elapsed, elapsed_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, p_err)
         build_samples(trial) = elapsed_global
      ENDDO
      build_time = median_of_three(build_samples)

      allocate(push)
      directory_temp_bytes = 0
      IF (use_directory) THEN
         CALL build_directory_single(nlocal, ids_me, nrequest, ids_req, push, directory_temp_bytes)
      ELSE
         CALL build_worker_pushdata(nlocal, ids_me, nrequest, ids_req, push)
      ENDIF
      CALL MPI_Allreduce(directory_temp_bytes, directory_temp_bytes_global, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, p_err)
      DO iteration = 1, 4
         CALL worker_push_data(push, send, recv, -1._r8)
      ENDDO

      DO trial = 1, timing_trials
         CALL MPI_Barrier(MPI_COMM_WORLD, p_err)
         t0 = MPI_Wtime()
         DO iteration = 1, push_repeats_used
            CALL worker_push_data(push, send, recv, -1._r8)
         ENDDO
         elapsed = (MPI_Wtime() - t0) / real(push_repeats_used, r8)
         CALL MPI_Allreduce(elapsed, elapsed_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, p_err)
         push_samples(trial) = elapsed_global
      ENDDO
      push_time = median_of_three(push_samples)

      failures = 0
      IF (nrequest > 0) failures = count(recv /= real(ids_req, r8))
      CALL MPI_Allreduce(failures, failures_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, p_err)
      benchmark_failures = benchmark_failures + failures_global

      IF (rank == 0) write(*,'(A,",",A,",",I0,",",I0,",",I0,",",I0,",",I0,",",ES14.6,",",ES14.6,",",I0,",",I0)') &
         trim(implementation_name(use_directory)), trim(label), nranks, base_nlocal, nrequest, &
         build_repeats_used, push_repeats_used, build_time, push_time, directory_temp_bytes_global, failures_global

      deallocate(push, ids_me, ids_req, send, recv)
   END SUBROUTINE run_single_case


   SUBROUTINE run_multi_case(label, base_nlocal, nrequest, build_repeats, push_repeats, empty_rank, use_directory)
      USE MOD_WorkerPushData, only: worker_pushdata_type, build_worker_pushdata, worker_push_data
      character(len=*), intent(in) :: label
      integer, intent(in) :: base_nlocal, nrequest, build_repeats, push_repeats
      logical, intent(in) :: empty_rank, use_directory

      type(worker_pushdata_type), allocatable :: push
      integer, allocatable :: ids_me(:), ids_req_single(:), ids_req(:,:)
      real(r8), allocatable :: area_req(:,:), send(:), recv(:), expected(:)
      integer :: nlocal, iteration, trial, failures, failures_global
      integer :: build_repeats_used, push_repeats_used
      integer :: directory_temp_bytes, directory_temp_bytes_global
      real(r8) :: t0, elapsed, elapsed_global
      real(r8) :: build_samples(timing_trials), push_samples(timing_trials)
      real(r8) :: build_time, push_time

      CALL build_case_ids(base_nlocal, nrequest, empty_rank, nlocal, ids_me, ids_req_single)
      allocate(ids_req(3, nrequest), area_req(3, nrequest), send(nlocal), recv(nrequest), expected(nrequest))
      IF (nlocal > 0) send = real(ids_me, r8)
      IF (nrequest > 0) THEN
         ids_req(1,:) = ids_req_single
         ids_req(2,:) = cshift(ids_req_single, 1)
         ids_req(3,:) = ids_req_single
         area_req(1,:) = 1._r8
         area_req(2,:) = 2._r8
         area_req(3,:) = 3._r8
         expected = 4._r8 * real(ids_req(1,:), r8) + 2._r8 * real(ids_req(2,:), r8)
      ENDIF
      CALL benchmark_repeat_counts(build_repeats, push_repeats, build_repeats_used, push_repeats_used)

      DO trial = 1, timing_trials
         CALL MPI_Barrier(MPI_COMM_WORLD, p_err)
         t0 = MPI_Wtime()
         DO iteration = 1, build_repeats_used
            allocate(push)
            IF (use_directory) THEN
               CALL build_directory_multi(nlocal, ids_me, nrequest, ids_req, area_req, push, directory_temp_bytes)
            ELSE
               CALL build_worker_pushdata(nlocal, ids_me, nrequest, ids_req, area_req, push)
            ENDIF
            deallocate(push)
         ENDDO
         elapsed = (MPI_Wtime() - t0) / real(build_repeats_used, r8)
         CALL MPI_Allreduce(elapsed, elapsed_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, p_err)
         build_samples(trial) = elapsed_global
      ENDDO
      build_time = median_of_three(build_samples)

      allocate(push)
      directory_temp_bytes = 0
      IF (use_directory) THEN
         CALL build_directory_multi(nlocal, ids_me, nrequest, ids_req, area_req, push, directory_temp_bytes)
      ELSE
         CALL build_worker_pushdata(nlocal, ids_me, nrequest, ids_req, area_req, push)
      ENDIF
      CALL MPI_Allreduce(directory_temp_bytes, directory_temp_bytes_global, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, p_err)
      DO iteration = 1, 4
         CALL worker_push_data(push, send, recv, -1._r8, mode='sum')
      ENDDO

      DO trial = 1, timing_trials
         CALL MPI_Barrier(MPI_COMM_WORLD, p_err)
         t0 = MPI_Wtime()
         DO iteration = 1, push_repeats_used
            CALL worker_push_data(push, send, recv, -1._r8, mode='sum')
         ENDDO
         elapsed = (MPI_Wtime() - t0) / real(push_repeats_used, r8)
         CALL MPI_Allreduce(elapsed, elapsed_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, p_err)
         push_samples(trial) = elapsed_global
      ENDDO
      push_time = median_of_three(push_samples)

      failures = 0
      IF (nrequest > 0) failures = count(abs(recv - expected) > 1.e-10_r8)
      CALL MPI_Allreduce(failures, failures_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, p_err)
      benchmark_failures = benchmark_failures + failures_global

      IF (rank == 0) write(*,'(A,",",A,",",I0,",",I0,",",I0,",",I0,",",I0,",",ES14.6,",",ES14.6,",",I0,",",I0)') &
         trim(implementation_name(use_directory)), trim(label), nranks, base_nlocal, nrequest, &
         build_repeats_used, push_repeats_used, build_time, push_time, directory_temp_bytes_global, failures_global

      deallocate(push, ids_me, ids_req_single, ids_req, area_req, send, recv, expected)
   END SUBROUTINE run_multi_case


   FUNCTION implementation_name(use_directory) RESULT(name)
      logical, intent(in) :: use_directory
      character(len=9) :: name

      IF (use_directory) THEN
         name = 'directory'
      ELSE
         name = 'ring'
      ENDIF
   END FUNCTION implementation_name


   SUBROUTINE build_directory_single(num_me, ids_me, num_req, ids_req, push, temp_bytes)
      USE MOD_WorkerPushData, only: worker_pushdata_type
      integer, intent(in) :: num_me, ids_me(:), num_req, ids_req(:)
      type(worker_pushdata_type), intent(inout) :: push
      integer, intent(out) :: temp_bytes

      integer :: i, iloc, n_req_uniq
      integer, allocatable :: ids_req_uniq(:)

      n_req_uniq = 0
      allocate(ids_req_uniq(num_req))
      DO i = 1, num_req
         CALL insert_into_sorted_list1(ids_req(i), n_req_uniq, ids_req_uniq, iloc)
      ENDDO

      IF (num_req > 0) THEN
         allocate(push%addr_single(num_req))
         DO i = 1, num_req
            push%addr_single(i) = find_in_sorted_list1(ids_req(i), n_req_uniq, ids_req_uniq(1:n_req_uniq))
         ENDDO
      ENDIF

      push%num_req_uniq = n_req_uniq
      push%mapping_kind = mapping_single
      CALL build_directory_uniq(num_me, ids_me, n_req_uniq, ids_req_uniq(1:n_req_uniq), push, temp_bytes)
      deallocate(ids_req_uniq)
   END SUBROUTINE build_directory_single


   SUBROUTINE build_directory_multi(num_me, ids_me, num_req, ids_req, area_req, push, temp_bytes)
      USE MOD_WorkerPushData, only: worker_pushdata_type
      integer, intent(in) :: num_me, ids_me(:), num_req, ids_req(:,:)
      real(r8), intent(in) :: area_req(:,:)
      type(worker_pushdata_type), intent(inout) :: push
      integer, intent(out) :: temp_bytes

      integer :: i, j, iloc, ndim1, n_req_uniq, owner
      integer, allocatable :: ids_req_uniq(:)
      logical, allocatable :: id_found(:)

      ndim1 = size(ids_req, 1)
      n_req_uniq = 0
      allocate(ids_req_uniq(ndim1 * num_req))
      DO j = 1, num_req
         DO i = 1, ndim1
            CALL insert_into_sorted_list1(ids_req(i,j), n_req_uniq, ids_req_uniq, iloc)
         ENDDO
      ENDDO

      IF (num_req > 0) THEN
         allocate(push%addr_multi(ndim1, num_req))
         DO j = 1, num_req
            DO i = 1, ndim1
               push%addr_multi(i,j) = &
                  find_in_sorted_list1(ids_req(i,j), n_req_uniq, ids_req_uniq(1:n_req_uniq))
            ENDDO
         ENDDO
      ENDIF

      push%num_req_uniq = n_req_uniq
      push%mapping_kind = mapping_multi
      CALL build_directory_uniq(num_me, ids_me, n_req_uniq, ids_req_uniq(1:n_req_uniq), push, temp_bytes)

      IF (num_req > 0) THEN
         allocate(push%area_multi(ndim1, num_req), push%sum_area(num_req))
         push%area_multi = area_req
         WHERE ((area_req <= 0._r8) .or. (ids_req <= 0))
            push%area_multi = 0._r8
         END WHERE

         allocate(id_found(n_req_uniq))
         id_found = .false.
         IF (push%nself > 0) id_found(push%self_to) = .true.
         DO owner = 0, nranks - 1
            IF (push%n_from_other(owner) > 0) id_found(push%other_to(owner)%val) = .true.
         ENDDO
         DO j = 1, num_req
            DO i = 1, ndim1
               IF (.not. id_found(push%addr_multi(i,j))) push%area_multi(i,j) = 0._r8
            ENDDO
         ENDDO
         push%sum_area = sum(push%area_multi, dim=1)
         deallocate(id_found)
      ENDIF

      deallocate(ids_req_uniq)
   END SUBROUTINE build_directory_multi


   SUBROUTINE build_directory_uniq(num_me, ids_me, n_req_uniq, ids_req_uniq, push, temp_bytes)
      USE MOD_WorkerPushData, only: worker_pushdata_type
      integer, intent(in) :: num_me, ids_me(:), n_req_uniq, ids_req_uniq(:)
      type(worker_pushdata_type), intent(inout) :: push
      integer, intent(out) :: temp_bytes

      integer :: owner, requester, i, iloc, istt, iend, idsp, total_owned
      integer :: total_send, total_recv, integer_bytes
      integer, allocatable :: owner_counts(:), owner_displs(:), directory_ids(:), directory_order(:)
      integer, allocatable :: source_loc(:,:), send_counts(:), recv_counts(:)
      integer, allocatable :: send_displs(:), recv_displs(:), send_locs(:), recv_locs(:)

      allocate(owner_counts(0:nranks-1), owner_displs(0:nranks-1))
      CALL MPI_Allgather(num_me, 1, MPI_INTEGER, owner_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, p_err)
      owner_displs(0) = 0
      DO owner = 1, nranks - 1
         owner_displs(owner) = owner_displs(owner-1) + owner_counts(owner-1)
      ENDDO
      total_owned = sum(owner_counts)

      allocate(directory_ids(total_owned), directory_order(total_owned))
      CALL MPI_Allgatherv(ids_me, num_me, MPI_INTEGER, directory_ids, owner_counts, owner_displs, &
         MPI_INTEGER, MPI_COMM_WORLD, p_err)

      ! Sort each owner's segment independently with the same helper used by
      ! production.  This preserves the one-match-per-owner behavior even if
      ! an ID is (unexpectedly) duplicated locally or across workers.
      DO owner = 0, nranks - 1
         IF (owner_counts(owner) > 0) THEN
            istt = owner_displs(owner) + 1
            iend = istt + owner_counts(owner) - 1
            directory_order(istt:iend) = [(i, i=1,owner_counts(owner))]
            CALL quicksort(owner_counts(owner), directory_ids(istt:iend), directory_order(istt:iend))
         ENDIF
      ENDDO

      allocate(source_loc(0:nranks-1, n_req_uniq))
      source_loc = 0
      DO owner = 0, nranks - 1
         IF (owner_counts(owner) <= 0) CYCLE
         istt = owner_displs(owner) + 1
         iend = istt + owner_counts(owner) - 1
         DO i = 1, n_req_uniq
            iloc = find_in_sorted_list1(ids_req_uniq(i), owner_counts(owner), directory_ids(istt:iend))
            IF (iloc > 0) source_loc(owner,i) = directory_order(istt+iloc-1)
         ENDDO
      ENDDO

      push%nself = count(source_loc(rank,:) > 0)
      IF (push%nself > 0) THEN
         allocate(push%self_from(push%nself), push%self_to(push%nself))
         push%self_from = pack(source_loc(rank,:), source_loc(rank,:) > 0)
         push%self_to = pack([(i, i=1,n_req_uniq)], source_loc(rank,:) > 0)
      ENDIF

      allocate(push%n_to_other(0:nranks-1), push%to_other(0:nranks-1))
      allocate(push%n_from_other(0:nranks-1), push%other_to(0:nranks-1))
      allocate(send_counts(0:nranks-1), recv_counts(0:nranks-1))
      DO owner = 0, nranks - 1
         send_counts(owner) = count(source_loc(owner,:) > 0)
      ENDDO
      send_counts(rank) = 0
      CALL MPI_Alltoall(send_counts, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, p_err)
      push%n_from_other = send_counts
      push%n_to_other = recv_counts

      DO owner = 0, nranks - 1
         IF (send_counts(owner) > 0) THEN
            allocate(push%other_to(owner)%val(send_counts(owner)))
            push%other_to(owner)%val = pack([(i, i=1,n_req_uniq)], source_loc(owner,:) > 0)
         ENDIF
      ENDDO

      allocate(send_displs(0:nranks-1), recv_displs(0:nranks-1))
      send_displs(0) = 0
      recv_displs(0) = 0
      DO owner = 1, nranks - 1
         send_displs(owner) = send_displs(owner-1) + send_counts(owner-1)
         recv_displs(owner) = recv_displs(owner-1) + recv_counts(owner-1)
      ENDDO
      total_send = sum(send_counts)
      total_recv = sum(recv_counts)
      allocate(send_locs(total_send), recv_locs(total_recv))
      DO owner = 0, nranks - 1
         idsp = send_displs(owner)
         DO i = 1, n_req_uniq
            IF (owner /= rank .and. source_loc(owner,i) > 0) THEN
               idsp = idsp + 1
               send_locs(idsp) = source_loc(owner,i)
            ENDIF
         ENDDO
      ENDDO
      CALL MPI_Alltoallv(send_locs, send_counts, send_displs, MPI_INTEGER, recv_locs, recv_counts, &
         recv_displs, MPI_INTEGER, MPI_COMM_WORLD, p_err)

      DO requester = 0, nranks - 1
         IF (recv_counts(requester) > 0) THEN
            allocate(push%to_other(requester)%val(recv_counts(requester)))
            istt = recv_displs(requester) + 1
            iend = istt + recv_counts(requester) - 1
            push%to_other(requester)%val = recv_locs(istt:iend)
         ENDIF
      ENDDO

      push%required_send_size = 0
      IF (push%nself > 0) push%required_send_size = maxval(push%self_from)
      DO requester = 0, nranks - 1
         IF (push%n_to_other(requester) > 0) &
            push%required_send_size = max(push%required_send_size, maxval(push%to_other(requester)%val))
      ENDDO

      integer_bytes = storage_size(0) / 8
      temp_bytes = integer_bytes * (2 * total_owned + nranks * n_req_uniq + 6 * nranks + total_send + total_recv)

      deallocate(owner_counts, owner_displs, directory_ids, directory_order, source_loc)
      deallocate(send_counts, recv_counts, send_displs, recv_displs, send_locs, recv_locs)
   END SUBROUTINE build_directory_uniq


   FUNCTION median_of_three(samples) RESULT(median)
      real(r8), intent(in) :: samples(timing_trials)
      real(r8) :: median

      median = sum(samples) - minval(samples) - maxval(samples)
   END FUNCTION median_of_three


   SUBROUTINE benchmark_repeat_counts(configured_build, configured_push, used_build, used_push)
      integer, intent(in) :: configured_build, configured_push
      integer, intent(out) :: used_build, used_push

      used_build = configured_build
      used_push = configured_push
      ! At 16+ ranks a single ring build already has enough duration for the
      ! timer.  Capping the inner loop keeps oversubscribed CI runs bounded;
      ! three independent trials still provide a robust median.
      IF (nranks >= 16) THEN
         used_build = min(used_build, 1)
         used_push = min(used_push, 4)
      ENDIF
   END SUBROUTINE benchmark_repeat_counts


   SUBROUTINE build_case_ids(base_nlocal, nrequest, empty_rank, nlocal, ids_me, ids_req)
      integer, intent(in) :: base_nlocal, nrequest
      logical, intent(in) :: empty_rank
      integer, intent(out) :: nlocal
      integer, allocatable, intent(out) :: ids_me(:), ids_req(:)
      integer :: i, owner, local_id, owner_count

      nlocal = base_nlocal
      IF (empty_rank .and. rank == 0) nlocal = 0
      allocate(ids_me(nlocal), ids_req(nrequest))

      IF (nlocal > 0) THEN
         IF (empty_rank) THEN
            ids_me = (rank - 1) * base_nlocal + [(i, i=1,nlocal)]
         ELSE
            ids_me = rank * base_nlocal + [(i, i=1,nlocal)]
         ENDIF
      ENDIF

      IF (nrequest <= 0) RETURN
      owner_count = nranks
      IF (empty_rank) owner_count = nranks - 1
      DO i = 1, nrequest
         IF (empty_rank) THEN
            owner = 1 + modulo(rank + i / 2, owner_count)
            local_id = 1 + modulo(17 * (i / 2) + rank, base_nlocal)
            ids_req(i) = (owner - 1) * base_nlocal + local_id
         ELSE
            owner = modulo(rank + i / 2, nranks)
            local_id = 1 + modulo(17 * (i / 2) + rank, base_nlocal)
            ids_req(i) = owner * base_nlocal + local_id
         ENDIF
      ENDDO
   END SUBROUTINE build_case_ids

END PROGRAM worker_push_mapping_benchmark
