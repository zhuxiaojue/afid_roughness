      module ibm_param
       implicit none
       integer :: npunx,npuny,npunz,npunte,mpun
       integer :: n
       parameter (mpun=1000000)
       real,allocatable,dimension(:,:,:) :: forclo
       real,allocatable,dimension(:,:) :: plth1, plth2
       integer,dimension(3,mpun,3) :: indgeo, indgeoe
       integer,dimension(mpun,3) :: indgeot, indgeoet
       real,dimension(3,mpun) :: distb
       real,dimension(mpun) :: distbt
       real,dimension(mpun) :: temb
       real,dimension(mpun) :: q1bo,q2bo,q3bo,densb

      end module ibm_param


