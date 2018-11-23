# ===========================================================================
#
# gds2bgen.r: format conversion between GDS and BGEN
#
# Copyright (C) 2018    Xiuwen Zheng (zhengxwen@gmail.com)
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# gds2bgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with gds2bgen.
# If not, see <http://www.gnu.org/licenses/>.


#############################################################
# Get bgen file information
#
seqBGEN_Info <- function(bgen.fn, verbose=TRUE)
{
    stopifnot(is.character(bgen.fn), length(bgen.fn)==1L)
    rv <- .Call(SEQ_BGEN_Info, bgen.fn)
    names(rv) <- c("num.sample", "num.variant", "compression", "layout",
        "sample.id")
    if (verbose)
    {
        cat("bgen file: ", normalizePath(bgen.fn), "\n", sep="")
        cat("# of samples: ", rv$num.sample, "\n", sep="")
        cat("# of variants: ", rv$num.variant, "\n", sep="")
        cat("compression method: ", rv$compression, "\n", sep="")
        cat("layout version: ", rv$layout, "\n", sep="")
       	cat("sample id: ")
        if (is.null(rv$sample.id))
        {
            cat("<anonymized>\n")
        } else {
            if (length(rv$sample.id) > 4L)
                cat(c(rv$sample.id[seq_len(4L)], "..."), sep=", ")
            else
                cat(rv$sample.id, sep=", ")
            cat("\n")
        }
    }
    invisible(rv)
}



#############################################################
# Format conversion from BGEN to GDS
#
seqBGEN2GDS <- function(bgen.fn, out.fn, storage.option="LZMA_RA",
    float.type=c("packed16", "packed8", "single", "double"),
    dosage=FALSE, prob=TRUE, start=1L, count=-1L,
    optimize=TRUE, digest=TRUE, parallel=FALSE, verbose=TRUE)
{
    stopifnot(is.character(bgen.fn), length(bgen.fn)==1L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    float.type <- match.arg(float.type)

    if (is.character(storage.option))
    {
        storage.option <- seqStorageOption(storage.option)
        s <- switch(float.type,
            packed16 = "packedreal16:offset=0,scale=0.0001",
            packed8  = "packedreal8:offset=0,scale=0.01",
            single   = "float32",
            double   = "float64")
		storage.option$mode <- c(
			`annotation/format/DS`=s, `annotation/format/GP`=s
		)
    }
    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))

    stopifnot(is.logical(dosage), length(dosage)==1L)
    stopifnot(is.logical(prob), length(prob)==1L)


    # get bgen info
    info <- seqBGEN_Info(bgen.fn, verbose=FALSE)

    if (verbose)
    {
        cat(date(), "\n", sep="")
        cat("BGEN Import:\n")
        cat("    file: ", bgen.fn, "\n", sep="")
        cat("    # of samples: ", info$num.sample, "\n", sep="")
        cat("    # of variants: ", info$num.variant, "\n", sep="")
        cat("    compression method: ", info$compression, "\n", sep="")
        cat("    layout version: ", info$layout, "\n", sep="")
       	cat("    sample id: ")
        if (is.null(info$sample.id))
        {
            cat("<anonymized>\n")
        } else {
            if (length(info$sample.id) > 4L)
                cat(c(info$sample.id[seq_len(4L)], "..."), sep=", ")
            else
                cat(info$sample.id, sep=", ")
            cat("\n")
        }
        flush.console()
    }

    nSamp <- info$num.sample
    if (nSamp <= 0L)
        stop("No sample in the bgen file.")

    # create a GDS file

    gfile <- createfn.gds(out.fn)
    on.exit({ if (!is.null(gfile)) closefn.gds(gfile) })
    if (verbose)
        cat("Output:\n    ", out.fn, "\n", sep="")

    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(gfile, "description")
    put.attr.gdsn(n, "bgen.version", info$layout)


    # add sample.id
    if (is.null(info$sample.id))
        info$sample.id <- seq_len(nSamp)
    SeqArray:::.AddVar(storage.option, gfile, "sample.id", info$sample.id,
        closezip=TRUE)

    # add variant.id
    SeqArray:::.AddVar(storage.option, gfile, "variant.id", storage="int32")

    # add position
    SeqArray:::.AddVar(storage.option, gfile, "position", storage="int32")

    # add chromosome
    SeqArray:::.AddVar(storage.option, gfile, "chromosome", storage="string")

    # add allele
    SeqArray:::.AddVar(storage.option, gfile, "allele", storage="string")

    # add a folder for genotypes
    varGeno <- addfolder.gdsn(gfile, "genotype")
    put.attr.gdsn(varGeno, "VariableName", "GT")
    put.attr.gdsn(varGeno, "Description", "Genotype")

    # add phase folder
    varPhase <- addfolder.gdsn(gfile, "phase")

    # add annotation folder
    varAnnot <- addfolder.gdsn(gfile, "annotation")

    # add id
    SeqArray:::.AddVar(storage.option, varAnnot, "id", storage="string")
    # add qual
    SeqArray:::.AddVar(storage.option, varAnnot, "qual", storage="float")
    # add filter
    varFilter <- SeqArray:::.AddVar(storage.option, varAnnot, "filter",
        storage="int32")

    # VCF INFO
    varInfo <- addfolder.gdsn(varAnnot, "info")

    # add the FORMAT field
    varFormat <- addfolder.gdsn(varAnnot, "format")
    if (dosage)
    {
        DS <- addfolder.gdsn(varFormat, "DS")
        put.attr.gdsn(DS, "Number", ".")
        put.attr.gdsn(DS, "Type", "Float")
        put.attr.gdsn(DS, "Description", "Estimated alternate allele dosage")
        SeqArray:::.AddVar(storage.option, DS, "data", storage="float",
            valdim=c(nSamp, 0L))
        SeqArray:::.AddVar(storage.option, DS, "@data", storage="int32",
            visible=FALSE)
    }
    if (prob)
    {
        GP <- addfolder.gdsn(varFormat, "GP")
        put.attr.gdsn(GP, "Number", ".")
        put.attr.gdsn(GP, "Type", "Float")
        put.attr.gdsn(GP, "Description", "Genotype probabilities")
        SeqArray:::.AddVar(storage.option, GP, "data", storage="float",
            valdim=c(nSamp, 0L))
        SeqArray:::.AddVar(storage.option, GP, "@data", storage="int32",
            visible=FALSE)
    }

    # call C function
    .Call(SEQ_BGEN_Import, bgen.fn, gfile$root, verbose)

    # add annotation folder
    addfolder.gdsn(gfile, "sample.annotation")

    # RLE-coded chromosome
    SeqArray:::.optim_chrom(gfile)
    # create hash
    SeqArray:::.DigestFile(gfile, digest, verbose)

    # close the GDS file
    closefn.gds(gfile)
    gfile <- NULL

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }
    if (optimize)
    {
        if (verbose)
        {
            cat("Optimize the access efficiency ...\n")
            flush.console()
        }
        cleanup.gds(out.fn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.fn))
}
