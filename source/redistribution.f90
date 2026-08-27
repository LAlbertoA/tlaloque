module density_redistribution_data
#ifdef MPIP
    implicit none

    integer, parameter :: max_send_blocks = 8
    integer, parameter :: meta_size = 14
    integer, parameter :: tag_base = 1000

    logical :: density_redistribution_ready = .false.

    integer :: nsend_rho_blocks = 0
    integer :: nrecv_rho_blocks = 0

    integer, allocatable :: send_dest(:)
    integer, allocatable :: send_nblock(:)
    integer, allocatable :: send_block_id(:)

    integer, allocatable :: send_src_i1(:), send_src_i2(:)
    integer, allocatable :: send_src_j1(:), send_src_j2(:)
    integer, allocatable :: send_src_k1(:), send_src_k2(:)

    integer, allocatable :: send_dst_i1(:), send_dst_i2(:)
    integer, allocatable :: send_dst_j1(:), send_dst_j2(:)
    integer, allocatable :: send_dst_k1(:), send_dst_k2(:)

    integer, allocatable :: recv_source(:)
    integer, allocatable :: recv_block_id(:)
    integer, allocatable :: recv_nblock(:)

    integer, allocatable :: recv_dst_i1(:), recv_dst_i2(:)
    integer, allocatable :: recv_dst_j1(:), recv_dst_j2(:)
    integer, allocatable :: recv_dst_k1(:), recv_dst_k2(:)

    type buffer_type
       real, allocatable :: data(:)
    end type buffer_type

    type(buffer_type), allocatable :: sendbuf(:)
    type(buffer_type), allocatable :: recvbuf(:)

#ifdef MPIP
    integer, allocatable :: send_req(:)
    integer, allocatable :: recv_req(:)
#endif

    integer, parameter :: max_phi_recv_blocks = 8
    integer, parameter :: phi_meta_size = 14
    integer, parameter :: phi_tag_base = 2000

    logical :: phi_redistribution_ready = .false.

    integer :: nsend_phi_blocks = 0
    integer :: nrecv_phi_blocks = 0

    integer, allocatable :: phi_send_dest(:)
    integer, allocatable :: phi_send_block_id(:)
    integer, allocatable :: phi_send_nblock(:)

    integer, allocatable :: phi_send_src_i1(:), phi_send_src_i2(:)
    integer, allocatable :: phi_send_src_j1(:), phi_send_src_j2(:)
    integer, allocatable :: phi_send_src_k1(:), phi_send_src_k2(:)

    integer, allocatable :: phi_send_dst_i1(:), phi_send_dst_i2(:)
    integer, allocatable :: phi_send_dst_j1(:), phi_send_dst_j2(:)
    integer, allocatable :: phi_send_dst_k1(:), phi_send_dst_k2(:)

    integer, allocatable :: phi_recv_source(:)
    integer, allocatable :: phi_recv_block_id(:)
    integer, allocatable :: phi_recv_nblock(:)

    integer, allocatable :: phi_recv_dst_i1(:), phi_recv_dst_i2(:)
    integer, allocatable :: phi_recv_dst_j1(:), phi_recv_dst_j2(:)
    integer, allocatable :: phi_recv_dst_k1(:), phi_recv_dst_k2(:)

	integer, allocatable :: phi_recv_src_i1(:), phi_recv_src_i2(:)
	integer, allocatable :: phi_recv_src_j1(:), phi_recv_src_j2(:)
	integer, allocatable :: phi_recv_src_k1(:), phi_recv_src_k2(:)

    type(buffer_type), allocatable :: phi_sendbuf(:)
    type(buffer_type), allocatable :: phi_recvbuf(:)

#ifdef MPIP
    integer, allocatable :: phi_send_req(:)
    integer, allocatable :: phi_recv_req(:)
#endif

#endif
end module density_redistribution_data