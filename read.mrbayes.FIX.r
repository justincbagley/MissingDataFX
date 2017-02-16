#!/usr/bin/env Rscript

################################### read.mrbayes.FIX.R ###################################

library(tools)
library(ape)		## This is an ips dependency, so load it first.
library(ips)		## ips is the main package we use in this script.

##--Create modified version of ips 'read.mrbayes' function, that fixes a bug in the code:
read.mrbayes.FIX <- function (file, digits = NULL){
    tr <- read.nexus(file)
    nt <- Ntip(tr)
    ni <- Nnode(tr)
    nn <- nt + ni
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    LEFT <- grep("\\[", X)
    tab <- extractMRBAYESstats(file, nn)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) {
        w <- LEFT == RIGHT
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1)) 
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    end <- semico[semico > i2][1]
    x <- X[(i2 + 1):end]
    x <- unlist(strsplit(x, "[,; \t]"))
    x <- x[nzchar(x)]
    TRANS <- matrix(x, ncol = 2, byrow = TRUE)
    TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
    n <- dim(TRANS)[1]
    start <- semico[semico > i2][1] + 1
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    tree <- gsub("^.*= *", "", tree)
    tree <- read.tree(text = tree)
    nt <- Ntip(tree)
    ni <- Nnode(tree)
    tree$edge.length <- 1:(Ntip(tree) + Nnode(tree) - 1)
    tree <- write.tree(tree, file = "")
    id <- unlist(strsplit(tree, ":"))[-1]
    id <- gsub("(^[[:digit:]]+)([[:punct:]].*$)", "\\1", id)
    id <- as.numeric(id)
    nn <- nt + ni
    tab <- tab[match(1:nn, c(id, nn))]
    tab <- t(tab)
    tr <- read.nexus(file)
    tr <- list(edge = tr$edge, edge.length = tr$edge.length, 
        Nnode = tr$Nnode, tip.label = tr$tip.label, mrbayes = tab)
    class(tr) <- ("phylo")
    attr(tr, "origin") <- file
    tr
}

