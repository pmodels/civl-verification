/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  Copyright (C) by Argonne National Laboratory.
 *     See COPYRIGHT in top-level directory.
 *
 */
/* Routine to schedule a ring exchange based allreduce. The algorithm is
 * based on Baidu's ring based allreduce. http://andrew.gibiansky.com/ */

#include "mpiimpl_cvl.h"

int MPIR_Allgatherv_intra_ring(const void *sendbuf,
                               MPI_Aint sendcount,
                               MPI_Datatype sendtype,
                               void *recvbuf,
                               const MPI_Aint * recvcounts,
                               const MPI_Aint * displs,
                               MPI_Datatype recvtype, MPIR_Comm * comm_ptr, MPIR_Errflag_t errflag);

int MPIR_Allreduce_intra_ring(const void *sendbuf, void *recvbuf, MPI_Aint count,
                              MPI_Datatype datatype, MPI_Op op,
                              MPIR_Comm * comm, MPIR_Errflag_t errflag)
{
    int mpi_errno = MPI_SUCCESS;
    int i, src, dst;
    int nranks, is_inplace, rank;
    MPI_Aint extent;
    MPI_Aint lb, true_extent;
    MPI_Aint *cnts, *displs;    /* Created for the allgatherv call */
    int send_rank, recv_rank, total_count;
    double *tmpbuf;
    double *drecvbuf = recvbuf;
    int tag;
    MPIR_Request *reqs[2];      /* one send and one recv per transfer */

    MPIR_CHKLMEM_DECL();

    is_inplace = (sendbuf == MPI_IN_PLACE);
    MPIR_COMM_RANK_SIZE(comm, rank, nranks);

    MPIR_Datatype_get_extent_macro(datatype, extent);
    MPIR_Type_get_true_extent_impl(datatype, &lb, &true_extent);
    extent = MPL_MAX(extent, true_extent);

    MPIR_CHKLMEM_MALLOC(cnts, nranks * sizeof(MPI_Aint));
    MPIR_CHKLMEM_MALLOC(displs, nranks * sizeof(MPI_Aint));

    for (i = 0; i < nranks; i++)
        cnts[i] = 0;

    total_count = 0;
    for (i = 0; i < nranks; i++) {
        cnts[i] = (count + nranks - 1) / nranks;
        if (total_count + cnts[i] > count) {
            cnts[i] = count - total_count;
            break;
        } else
            total_count += cnts[i];
    }

    displs[0] = 0;
    for (i = 1; i < nranks; i++)
        displs[i] = displs[i - 1] + cnts[i - 1];

    /* Phase 1: copy to tmp buf */
    if (!is_inplace) {
        mpi_errno = MPIR_Localcopy(sendbuf, count, datatype, drecvbuf, count, datatype);
        MPIR_ERR_CHECK(mpi_errno);
    }

    /* Phase 2: Ring based send recv reduce scatter */
    /* Need only 2 spaces for current and previous reduce_id(s) */
    MPIR_CHKLMEM_MALLOC(tmpbuf, count * extent);

    src = (nranks + rank - 1) % nranks;
    dst = (rank + 1) % nranks;

    for (i = 0; i < nranks - 1; i++) {
        recv_rank = (nranks + rank - 2 - i) % nranks;
        send_rank = (nranks + rank - 1 - i) % nranks;

        /* KJR: use a single MPI_Sendrecv to bypass unsupported nonblocking ops */
        /* mpi_errno = MPIC_Irecv(tmpbuf, cnts[recv_rank], datatype, src, MPIR_ALLREDUCE_TAG, comm, &reqs[0]); */
        /* MPIR_ERR_CHECK(mpi_errno); */

        /* mpi_errno = MPIC_Isend((char *) recvbuf + displs[send_rank] * extent, cnts[send_rank], */
        /*                        datatype, dst, MPIR_ALLREDUCE_TAG, comm, &reqs[1], errflag); */
        /* MPIR_ERR_CHECK(mpi_errno); */

        /* mpi_errno = MPIC_Waitall(2, reqs, MPI_STATUSES_IGNORE); */
        /* MPIR_ERR_CHECK(mpi_errno); */
        MPIC_Sendrecv(drecvbuf + displs[send_rank], cnts[send_rank], datatype,
                      dst, MPIR_ALLREDUCE_TAG, tmpbuf, cnts[recv_rank], datatype, src,
                      MPIR_ALLREDUCE_TAG, comm, MPI_STATUS_IGNORE, errflag);

        mpi_errno =
            MPIR_Reduce_local(tmpbuf, drecvbuf + displs[recv_rank],
                              cnts[recv_rank], datatype, op);
        MPIR_ERR_CHECK(mpi_errno);
    }

    /* Phase 3: Allgatherv ring, so everyone has the reduced data */
    mpi_errno = MPIR_Allgatherv_intra_ring(MPI_IN_PLACE, -1, MPI_DATATYPE_NULL, drecvbuf, cnts,
                                           displs, datatype, comm, errflag);
    MPIR_ERR_CHECK(mpi_errno);

  fn_exit:
    MPIR_CHKLMEM_FREEALL();
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Allgatherv_intra_ring(const void *sendbuf,
                               MPI_Aint sendcount,
                               MPI_Datatype sendtype,
                               void *recvbuf,
                               const MPI_Aint * recvcounts,
                               const MPI_Aint * displs,
                               MPI_Datatype recvtype, MPIR_Comm * comm_ptr, MPIR_Errflag_t errflag)
{
    int comm_size, rank, i, left, right;
    int mpi_errno = MPI_SUCCESS;
    MPI_Status status;
    MPI_Aint recvtype_extent;
    MPI_Aint total_count;
    double *drecvbuf = recvbuf;

    MPIR_COMM_RANK_SIZE(comm_ptr, rank, comm_size);

    total_count = 0;
    for (i = 0; i < comm_size; i++)
        total_count += recvcounts[i];

    if (total_count == 0)
        goto fn_exit;

    MPIR_Datatype_get_extent_macro(recvtype, recvtype_extent);

    if (sendbuf != MPI_IN_PLACE) {
        /* First, load the "local" version in the recvbuf. */
        mpi_errno = MPIR_Localcopy(sendbuf, sendcount, sendtype,
                                   recvbuf + displs[rank],
                                   recvcounts[rank], recvtype);
        MPIR_ERR_CHECK(mpi_errno);
    }

    left = (comm_size + rank - 1) % comm_size;
    right = (rank + 1) % comm_size;

    MPI_Aint torecv, tosend, max, chunk_count;
    torecv = total_count - recvcounts[rank];
    tosend = total_count - recvcounts[right];

    chunk_count = 0;
    max = recvcounts[0];
    for (i = 1; i < comm_size; i++)
        if (max < recvcounts[i])
            max = recvcounts[i];
    if (MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE > 0 &&
        max * recvtype_extent > MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE) {
        chunk_count = MPIR_CVAR_ALLGATHERV_PIPELINE_MSG_SIZE / recvtype_extent;
        /* Handle the case where the datatype extent is larger than
         * the pipeline size. */
        if (!chunk_count)
            chunk_count = 1;
    }
    /* pipeline is disabled */
    if (!chunk_count)
        chunk_count = max;

    MPI_Aint soffset, roffset;
    int sidx, ridx;
    sidx = rank;
    ridx = left;
    soffset = 0;
    roffset = 0;
    while (tosend || torecv) {  /* While we have data to send or receive */
        MPI_Aint sendnow, recvnow;
        sendnow = ((recvcounts[sidx] - soffset) >
                   chunk_count) ? chunk_count : (recvcounts[sidx] - soffset);
        recvnow = ((recvcounts[ridx] - roffset) >
                   chunk_count) ? chunk_count : (recvcounts[ridx] - roffset);

        double *sbuf, *rbuf;
        sbuf = recvbuf + (displs[sidx] + soffset);
        rbuf = recvbuf + (displs[ridx] + roffset);

        /* Protect against wrap-around of indices */
        if (!tosend)
            sendnow = 0;
        if (!torecv)
            recvnow = 0;

        /* Communicate */
        if (!sendnow && !recvnow) {
            /* Don't do anything. This case is possible if two
             * consecutive processes contribute 0 bytes each. */
        } else if (!sendnow) {  /* If there's no data to send, just do a recv call */
            mpi_errno =
                MPIC_Recv(rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG, comm_ptr, &status);
            MPIR_ERR_CHECK(mpi_errno);
            torecv -= recvnow;
        } else if (!recvnow) {  /* If there's no data to receive, just do a send call */
            mpi_errno =
                MPIC_Send(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG, comm_ptr, errflag);
            MPIR_ERR_CHECK(mpi_errno);
            tosend -= sendnow;
        } else {        /* There's data to be sent and received */
            mpi_errno = MPIC_Sendrecv(sbuf, sendnow, recvtype, right, MPIR_ALLGATHERV_TAG,
                                      rbuf, recvnow, recvtype, left, MPIR_ALLGATHERV_TAG,
                                      comm_ptr, &status, errflag);
            MPIR_ERR_CHECK(mpi_errno);
            tosend -= sendnow;
            torecv -= recvnow;
        }

        soffset += sendnow;
        roffset += recvnow;
        if (soffset == recvcounts[sidx]) {
            soffset = 0;
            sidx = (sidx + comm_size - 1) % comm_size;
        }
        if (roffset == recvcounts[ridx]) {
            roffset = 0;
            ridx = (ridx + comm_size - 1) % comm_size;
        }
    }

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}
