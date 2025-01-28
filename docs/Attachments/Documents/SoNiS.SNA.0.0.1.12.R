# SoNiS.SNA functions version 0.0.1.12
# 28.01.2025 
# These functions are part of a self-written package that I created during my masters thesis and
# developed through my PhD. I plan to release a full package at some point, but for now,
# the functions are presented and documented (at least to some extent) here.


#### Attribute to Matrix ####
# this function takes a data frame with 2 columns
# the first column is expected to be the identifiers for the actors,
# the second column is expected to be the attribute expression for the actors.
# By choosing a mode, the function turns the attribute data frame into a matrix format
# that can be used e.g. for QAP correlations or regressions
attribute_to_matrix <- function(df,mode="match"){
  stopifnot("Mode invalid or unknown"=mode%in%c("match","diff","abs_diff","sq_diff","prod","sum","receiver","sender"))
  stopifnot("Input is not a data frame or tibble"=is.data.frame(df))
  n <- length(df[[1]])
  matrix(0,nrow=n,ncol=n) %>%
    `dimnames<-`(list(df[[1]],df[[1]])) -> attr_mat
  if(mode=="receiver"){
    matrix(rep(df[[2]],times=n),nrow=n,byrow=TRUE) %>%
      `dimnames<-`(list(df[[1]],df[[1]])) -> attr_mat
  }
  if(mode=="sender"){
    matrix(rep(df[[2]],times=n),nrow=n,byrow=FALSE) %>%
      `dimnames<-`(list(df[[1]],df[[1]])) -> attr_mat
  }
  if(mode=="match"){
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if(df[i,2]==df[j,2]){
          attr_mat[[i,j]] <- 1
          attr_mat[[j,i]] <- 1
        }
      }
    }
  }else{
    for(i in 1:n){
      for(j in 1:n){
        if(mode=="diff"){
          attr_mat[[i,j]] <- (df[i,2]-df[j,2])
        }
        if(mode=="abs_diff"){
          attr_mat[[i,j]] <- abs(df[i,2]-df[j,2])
        }
        if(mode=="sq_diff"){
          attr_mat[[i,j]] <- (df[i,2]-df[j,2])^2
        }
        if(mode=="prod"){
          attr_mat[[i,j]] <- (df[i,2]*df[j,2])
        }
        if(mode=="sum"){
          attr_mat[[i,j]] <- (df[i,2]+df[j,2])
        }
      }
    }
  }
  attr_mat
}

#### Data to Adjacency ####
# convert directly from data to adjacency matrix / network object
# this function essentially chains data_to_nodelist and nodelist_to_adjacency
data_to_adjacency <- function(data,school,class,item,participants="all",as_network=TRUE,...){
  stopifnot("Choice of participants not available or outdated" = participants %in% c("only","all"))
  # Filter and select data to display the nodelist for the selected school and class
  nodelist <- data_to_nodelist(data,school,class,item)
  # Convert nodelist to adjacency matrix
  adj <- nodelist_to_adjacency(data = nodelist,
                                          as_network = as_network,
                                          participants = participants)
  
  adj
}

#### Data to Nodelist ####
# this is essentially a filtering and selection function that mainly works with our datasets
# because every class has a school (SC) and class (KC) identifier.
# this function is only usable if you also use these columns to identify your classes / networks
# but is is essentially just a nice wrapper to filter the data
data_to_nodelist <- function(data, school, class, item){
  stopifnot("Data is not a data.frame, tibble, matrix or list"= typeof(data) %in% c("data.frame","tibble","matrix","list"))
  if(typeof(data)=="matrix"){
    data <- data.frame(data)
  }
  data %>%
    dplyr::filter(SC==school , KC==class) %>%
    dplyr::select(c("ZC",dplyr::contains(item))) -> res
  res
}

#### Selection Index ####
# this function chooses which actors to include in the adjacency matrix
# "only" includes only the actors that were present in the data collection 
# (meaning they have a row in the dataset)
# "all" includes also those who didn't participate but were nominated by others
# these students then don't have outgoing ties but the incoming ties are counted
selection_index <- function(data,participants){
  # Vorher stoppen wenn falsche Parameter eingegeben wurden
  stopifnot("Datentyp unbekannt oder unkompatibel." = typeof(data) %in% c("data.frame","tibble","list","matrix","integer","double"))
  stopifnot("Participants Auswahlmethode unbekannt" = participants %in% c("only","all"))
  # Datentyp zu Matrix konvertieren
  if(typeof(data)=="data.frame" | typeof(data)=="tibble" | typeof(data)=="list"){
    mat <- as.matrix(data)
  }else{
    mat <- data
  }
  if(participants=="only"){index_used <- sort(mat[,1])}
  if(participants=="all"){index_used <- sort(unique(as.vector(mat)))}
  index_used
}

#### Jaccard ####
# This function calculates the jaccard index between 2 networks or individuals in 2 networks
jaccard <- function(nets,net2=NULL,mode="index"){
  modes <- c("index","distance","individual_index_out","individual_distance_out","individual_index_in","individual_distance_in")
  stopifnot("Mode not known"=mode %in% modes)
  if(is.null(net2)){
    if(typeof(nets)=="list"){
      net1 <- nets[[1]]
      net2 <- nets[[2]]
    }
    if(is.array(nets)){
      net1 <- nets[1,,]
      net2 <- nets[2,,]
    }
  }else{
    net1 <- nets
  }
  net1 <- as.matrix(net1)
  net2 <- as.matrix(net2)
  stopifnot("Matrices are not the same size"=ncol(net1)==ncol(net2))
  if(mode=="index"){
    A <- c(net1)
    B <- c(net2)
    AB <- A+B
    C <- length(AB[which(AB==2)])
    N1 <- length(A[which(A==1)])
    N2 <- length(B[which(B==1)])
    res <- C/(N1+N2-C)
  }
  if(mode=="distance"){
    res <- (1-jaccard(net1,net2,mode="index"))
  }
  n <- nrow(net1)
  names <- dimnames(net1)[[1]]
  if(mode=="individual_index_out"){
    temp <- vector(mode="numeric",length(n))
    for(i in 1:n){
      temp[i] <- jaccard(net1[i,],net2[i,],mode="index")
    }
    names(temp) <- names
    res <- temp
  }
  if(mode=="individual_distance_out"){
    temp <- vector(mode="numeric",length(n))
    for(i in 1:n){
      temp[i] <- (1-jaccard(net1[i,],net2[i,],mode="index"))
    }
    names(temp) <- names
    res <- temp
  }
  if(mode=="individual_index_in"){
    temp <- vector(mode="numeric",length(n))
    for(i in 1:n){
      temp[i] <- jaccard(net1[,i],net2[,i],mode="index")
    }
    names(temp) <- names
    res <- temp
  }
  if(mode=="individual_distance_in"){
    temp <- vector(mode="numeric",length(n))
    for(i in 1:n){
      temp[i] <- (1-jaccard(net1[,i],net2[,i],mode="index"))
    }
    names(temp) <- names
    res <- temp
  }
  res
}

#### Multinet ####
# this function takes a list of networks and combines them to one matrix in a block-diagonal structure
# there is also the option to standardize the network before or after combining them
multinet <- function(netlist,scale=c(TRUE,FALSE),mode=c("before","after")){
  netlist <- lapply(netlist,as.matrix)
  numnets <- length(netlist)
  nvec <- vector(mode="numeric")
  for(i in 1:numnets){
    nvec[[i]] <- nrow(as.matrix(netlist[[i]]))
  }
  if(isTRUE(scale)){if(mode=="before"){netlist <- lapply(netlist,scale_net)}}
  meganet <- matrix(NA,nrow=sum(nvec),ncol=sum(nvec))
  meganet[(1:nvec[1]),(1:nvec[1])] <- netlist[[1]]
  for(j in 2:numnets){
    meganet[(sum(nvec[1:(j-1)])+1):(sum(nvec[1:j])),(sum(nvec[1:(j-1)])+1):(sum(nvec[1:j]))] <- netlist[[j]]
  }
  if(isTRUE(scale)){if(mode=="after"){meganet <- scale_net(meganet)}}
  meganet
}

#### Netlm Scaled ####
# this function is a wrapper for sna::netlm() to simplify the process to calculate standardized beta coefficients
# it takes 1 network as dependent variable (AV) 
# and a list of predictors (UV)
# you can also pass some additional arguments to netlm via "..." if you want to
# everything is scaled before the netlm resulting in standardized beta as a result of the regression
netlm_scaled <- function(AV,UV,...){
  stopifnot("AV must be a matrix or network"= is.matrix(AV) | class(AV)=="network")
  stopifnot("UV must be matrix, network or list of matrices or networks"= is.matrix(UV) | class(UV)=="network" | typeof(UV)=="list")
  AV <- scale_net(AV)
  if(class(UV)=="list"){
    UV <- lapply(UV,scale_net)
    
  }else{
    UV <- scale_net(UV)
  }
  res <- sna::netlm(AV,UV,...)
  res
}

#### Ngraph ####
# this is a simple wrapper for sna::gplot()
# it has some presets for visualization that I like
# it also has the option to color in reciprocal edges via the "reciprocal" and color options 
ngraph <- function(data,
                   label = network::network.vertex.names(data),
                   label.pos = 5,
                   vertex.col = "dark red",
                   label.col = "white",
                   vertex.cex = 1.2,
                   arrowhead.cex = 0.5,
                   reciprocal=FALSE,
                   reciprocal_color="red",
                   non_reciprocal_color="black",
                   ...){
  if(reciprocal==FALSE){
    sna::gplot(dat = data,
               label = label,
               label.pos = label.pos,
               vertex.col = vertex.col,
               label.col = label.col,
               vertex.cex = vertex.cex,
               arrowhead.cex = arrowhead.cex,
               ...)
  }
  if(reciprocal==TRUE){
    data_edgelist <- as.edgelist(data)
    n <- nrow(data_edgelist)
    reciprocal_edges <- vector(mode="numeric",length=n)
    for(i in 1:n){
      for(j in 1:n){
        if(all(data_edgelist[i,]==rev(data_edgelist[j,]))){
          reciprocal_edges[[i]] <- 1
          reciprocal_edges[[j]] <- 1
        }else{
          next
        }
      }
    }
    reciprocal_edges_attr <- ifelse(reciprocal_edges==0,non_reciprocal_color,reciprocal_color)
    sna::gplot(dat = data_edgelist,
               label = label,
               label.pos = label.pos,
               vertex.col = vertex.col,
               label.col = label.col,
               vertex.cex = vertex.cex,
               arrowhead.cex = arrowhead.cex,
               edge.col=reciprocal_edges_attr,
               ...)
  }
}

#### Nodelist to Adjacency ####
# this function takes a nodelist and converts it to an adjacency matrix or a network object
nodelist_to_adjacency <- function(data,as_network=TRUE,participants="all", netsize = 50, ...){
  # test for errors
  stopifnot("Datentyp unbekannt oder unkompatibel." = typeof(data)%in%c("data.frame","tibble","matrix","list"))
  stopifnot("Participants Auswahlmethode unbekannt oder veraltet." = participants %in% c("only","all"))
  # convert nodelist to matrix format
  if(typeof(data)=="data.frame" | typeof(data)=="tibble" | typeof(data)=="list"){
    mat <- as.matrix(data)
  }else{
    mat <- data
  }
  # choose index vector based on participant option
  index_used <- selection_index(data=mat,participants=participants)
  stopifnot("Nur eine Person ausgewählt!"=length(index_used)>1)
  # initialize matrix of the right size (standard is 50 actors)
  adj_start <- matrix(rep(0,times = netsize * netsize),nrow=netsize)
  # write in 1 for every nomination that is in the nodelist
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      adj_start[mat[i,1],mat[i,j]] <- 1
    }
  }
  # delete all rows that were not used (because in our dataset not every number is taken)
  adj_clean_rows <- adj_start[index_used,]
  adj_clean_all <- adj_clean_rows[,index_used]
  # output network object or adjacency matrix
  if(as_network==TRUE){
    res <- network::network(adj_clean_all,matrix.type="adjacency",...)
    network::network.vertex.names(res) <- index_used
  }else{
    res <- adj_clean_all
    dimnames(res) <- list(index_used,index_used)
  }
  res
}

#### Scale Net ####
# this function scales a network to mean = 0 and sd = 1
# it essentially deletes the diagonal of the matrix (because we exclude loops)
# then scales the rest of the network and reinstates the 0 on the diagonal
scale_net <- function(net){
  stopifnot("Object must be a matrix or network"= is.matrix(net) | class(net)=="network")
  netmat <- as.matrix(net)
  stopifnot("Network must not be empty"=sum(netmat,na.rm=TRUE)!=0)
  for(i in 1:nrow(netmat)){netmat[i,i] <- NA}
  scaled <- matrix(scale(c(netmat)),nrow=nrow(netmat))
  for(i in 1:nrow(scaled)){scaled[i,i] <- 0}
  scaled
}

#### Wahl und Ablehnung ####
# this function displays a couple of plots and metrics based on Petillon (1978) "Der unbeliebte Schüler"
# you need a positive (e.g., like) and a negative (e.g., dislike) network
# it then calculates sociometric positions for every actor and plots them
# as the original text is in German the comments here are also in German
wahl_und_ablehnung <- function(network_pos,network_neg,xlim = c(0.65,1.6), ylim = c(0.65, 1.6),color="black",base_size = 11,text_size=3,labels=TRUE){
  # Indegrees positiv und negativ bestimmen
  indeg_pos <- sna::degree(network_pos,cmode="indegree")
  indeg_neg <- sna::degree(network_neg,cmode="indegree")
  # Indegrees und ZC zuordnen
  names(indeg_pos) <- network::network.vertex.names(network_pos)
  names(indeg_neg) <- network::network.vertex.names(network_neg)
  # Wahl- und Ablehnungsstatus ausrechnen
  WS <- (1 + (indeg_pos-mean(indeg_pos))/(length(indeg_pos)-1))
  AS <- (1 + (indeg_neg-mean(indeg_neg))/(length(indeg_neg)-1))
  # Verteilung des Wahl- und Ablehnungsstatusses plotten
  # TODO: GGplot
  # Wahl- und Ablehnungsstatus für Einteilung in einen Dataframe zusammenfassen
  WS_AS <- tibble::tibble(ZC=names(WS),AS=AS,WS=WS)
  # Linien und Bezeichnungen der Typen nach Petillon 1980 vorbereiten
  type_lines <- tibble::tibble(X=c(0,2,0.8,0.8,0.8,0,2,1.2,1.2,1.2),
                               Y=c(0.8,0.8,0.8,2,1.2,1.2,1.2,1.2,2,0))
  type_labels <- tibble::tibble(X=c(1.4,1.4,1,1,0.68,1,0.68),
                                Y=c(0.68,1,0.68,1,1,1.4,1.4),
                                label=c("Typ 1: Ausgestoßener","Typ 2: Abgelehnter","Typ 3: Unbeachteter","Typ 4: Unauffälliger","Typ 5: \n Anerkannter","Typ 6: Beachteter","Typ 7: \n Star"))
  if(labels==TRUE){
    baseplot <- ggplot2::ggplot(data = WS_AS,mapping = ggplot2::aes(x = AS, y = WS, label = ZC))+
      ggrepel::geom_text_repel()
  }
  if(labels==FALSE){baseplot <- ggplot2::ggplot(data = WS_AS,mapping = ggplot2::aes(x = AS, y = WS))}
  plot <- baseplot +
    ggplot2::geom_point(color=color,show.legend = TRUE)+
    ggplot2::geom_path(data = type_lines, mapping = ggplot2::aes(x = X, y = Y, label = NULL),color="#555555")+
    ggplot2::geom_text(data = type_labels, mapping = ggplot2::aes(x = X, y = Y, label = label),color="#555555",size=text_size)+
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)+
    ggplot2::theme_bw(base_size = base_size)
  res <- list(indeg_pos=indeg_pos,
              indeg_neg=indeg_neg,
              WS=WS,
              AS=AS,
              WS_AS=WS_AS,
              plot=plot)
  res
}






