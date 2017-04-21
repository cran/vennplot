#' Draw Venn diagram in 2D or 3D
#'
#' @param combinations Named numeric vector or data.frame where each row is a logical vector indicating set membership.  See Details.
#' @param fit Scaling of semi-diameters; a scalar between 0.5 and 1.
#' @param GRAM If \code{TRUE}, use the centralized GRAM matrix as the initial location.  Otherwise use the center of an arbitrary circle.
#' @param arbitrary When \code{GRAM = FALSE}, an integer designating which circle to use for inital location, or \code{NULL} to use a random circle.  Argument is ignored when \code{GRAM = TRUE}.
#' @param main Title of 2D plot.
#' @param NAME The name of each circle. See Examples.
#' @param ALPHA Size of steepest descent.
#' @param alpha.col Color darkness.
#' @param ThreeD Draw Venn diagram in 3D. See Examples.
#' @param delta Closeness betwen connected components.
#' @param mu Generates defaults for unspecificed two-way intersections.  See Details.
#' @param disjoint If \code{TRUE} plots only the disjoint part of overlapping circles.
#' @param mar Plot margins.
#' @param priority Scalar or vector specifying the fitting priority for intersections.  See Examples.
#' @param weight Priority weights.  Must the the same length as \code{priority}.  See Examples.
#' @details The names of data sets (circles) supplied to \code{combination} must be single letters, e.g., \code{combination = c(A=1, B=2, AB=0.5)}.  Here \code{A} means the whole data-set, so does \code{B}, which means \code{AB} can be no larger than \code{min(A,B)}.
#' If a few two way intersections are unspecified, \code{mu} gives one way to generate them. For example, Suppose input combinations are \code{c(a=1,b=2,c=1,abc = 0.2)}.  Then default values for the two-way interserctions are  \code{ab = mu^(3-2)*abc = 0.2mu^(3-2)}, and \code{bc = mu^(3-2)*abc = 0.2mu^(3-2), ac = mu^(3-2)*abc = 0.2mu^(3-2)}.
#'
#'
#' @return An object of the class \code{vennplot} with following components:
#' \describe{
#'   \item{center}{centers of the circles (columns are (\code{x}, \code{y}) or (\code{x}, \code{y}, \code{z}) coordinates).}
#'   \item{semidiameters}{semi-diameters of the circles.}
#'   \item{LOSS}{total loss of \code{vennplot}.}
#'   \item{weighted.least.square}{Given specific priorities and weights, the weighted least square between input and lay-out.}
#' }
#' @examples
#' # arbitray sets
#' combinations = c(A=1.8, B=0.9,C=1.3, D = 1.3,E = 1.6, AC=0.3,
#'                  AD= 0.3,BE = 0.3, AE = 0.4, f = 0.7,g =0.8,
#'                  h = 0.5, gf = 0.2, Bh = 0.1,i = 1,j = 0.4,
#'                  k=0.7,l = 1.4,kl = 0.2,m = 0.5,lm = 0.2,
#'                  o = 0.8,p = 0.9, op = 0.3)
#' ve = vennplot(combinations)
#'
#' # named sets
#' # combinations = c(A=1, B=1,C=1, ABC = 0.1)
#' # ve = vennplot(combinations,NAME = c("Loon","Goose","Duck"))
#'
#' # effect of parameter mu
#' #combinations = c(A=1, B=1, C=1, D=1, E=1, ABCDE = 0.1)
#' # par(mfrow = c(1,2))
#' # ve = vennplot(combinations)
#' # ve = vennplot(combinations,mu=1.2)
#'
#' # 3D Venn plot
#' combinations = c(A=803, B=304,C=1015, D = 1100,E = 1005,f = 967,H=3020,
#'                  CD = 1000,BC = 248,ABC = 185,ADE = 327,CDfH = 846,
#'                  I=800, J=760,K=1000, L = 1100,M = 900, IK=333,
#'                  JM = 251,IL= 289, KM = 412,JL = 213)
#' timestart <- Sys.time()
#' ve = vennplot(combinations, disjoint = TRUE,ThreeD = TRUE)
#' timeend <- Sys.time()
#' runningtime <- timeend-timestart
#' print(runningtime)
#'
#' # effect of parameters priority and weight
#' # combinations = c(A = 79,B=29,C=58,AB = 25,AC = 44, BC=18,ABC = 1)
#' # par(mfrow = c(2,2))
#' # vennplot(combinations, priority = 2)
#' # vennplot(combinations, priority = 3)
#' # vennplot(combinations, priority = c(2,3), weight = c(30,1000))
#' # vennplot(combinations, priority = c(2,3), weight = c(30,10))
#'
#' # binary data
#' # combinations = sharks[,c(1,3:5,8)]
#' # vennplot(combinations = combinations)
#' @export
vennplot = function(combinations = NULL,fit = 0.5, GRAM = T,main = NULL, arbitrary = NULL, mu = 2,
                    NAME = NULL,ALPHA = 10^(-5),alpha.col = 0.3, disjoint = F,ThreeD = FALSE,
                    delta = 0.01, mar = rep(1,4),priority = 2,weight = rep(1,length(priority))){
  if(is.null(combinations)){
    stop("combinations should not be empty")
  }else if (is.list(combinations) || is.null(dim(combinations)) == F){
    if(is.null(NAME)){
      NAME = datatocom(combinations)$Name
    }
    combinations = datatocom(combinations)$combinations
  }
  if(disjoint == TRUE){
    combinations = disjoint.transform(combinations)
  }
  name = names(combinations)
  mc = str_length(name)
  if(max(mc)==2){
    combinations = combinations[which(combinations!=0)]
  }
  name = names(combinations)
  mc = str_length(name)
  m = length(which(mc==1))
  cols <- rainbow(m,alpha = alpha.col)
  name.dis = name[which(mc==1)]
  name.joint = name[-which(mc==1)]
  if(m==1){
    radius = 0.5
    lossnum = 0
    weighted.least.square = 0
    if(ThreeD){
      xy = matrix(rep(0,3),nrow = 1)
      colnames(xy) = c("x","y","z")
      if (is.null(NAME)){name = name.dis}else{name = NAME}
      open3d()
      spheres3d(xy[,1],xy[,2],xy[,3], radius = radius,
                color = cols, alpha = alpha.col)
      text3d(xy[,1],xy[,2],xy[,3], name)
    }
    else{
      xy = matrix(rep(0,2),nrow = 1)
      colnames(xy) = c("x","y")
      if (is.null(NAME)){name = name.dis}else{name = NAME}
      plot.new()
      circle.plot(xy,radius,name = name, col = cols, main = main, mar = mar)
    }
    outfinal = list(center = xy,diameter = radius, LOSS = lossnum, weighted.least.square = weighted.least.square)
  }
  else if(length(name.joint)==0){
    outev3 = EV3(combinations = combinations, m=m, ThreeD = ThreeD, cols = cols, main = main, mar = mar,
                 fit = fit, name.dis = name.dis,delta = delta, NAME=NAME, alpha.col = alpha.col)
    xy = outev3$xy
    radius = outev3$radius
    lossnum = 0
    weighted.least.square = 0
    outfinal = list(center = xy,diameter = radius, LOSS = lossnum, weighted.least.square = weighted.least.square)
  }
  else{
    combinations = combinations[order(mc)]
    name = names(combinations)
    mc = sort(mc)
    joint = combinations[-which(mc==1)]
    disjointcom = combinations[which(mc==1)]
    if(length(weight)!=length(priority)){stop("length of priority and length of weight should be the same")}
    if(ThreeD){priority = 2}
    if(length(priority)>1 && any(priority==1)){weight = weight[-which(priority==1)];priority = priority[which(priority!=1)]
    }else if(length(priority)==1 && priority==1){priority = 2}
    if(all(priority%in%mc)==FALSE){weight = weight[-which(priority%in%mc!=TRUE)];priority = priority[which(priority%in%mc == TRUE)]}
    if(length(priority)==0){priority = 2}
    if(length(priority)>1){mu = 1.1
    }else if(length(priority)==1 && priority!=2){mu = 1.1}
    if(max(mc)>=3){
      out = MC(combinations = combinations,disjointcom = disjointcom,mu = mu)
      if(length(priority)==1 && priority==2){
        name = out$name
        name.joint = out$name.joint
        joint = out$joint
        newcombinations = out$combinations
        newmc = out$mc
      }else{
        out1 = STC(mu = mu, newcombinations = out$combinations, priority = priority,disjointcom = disjointcom,
                   combinations = combinations, group = NULL, plus = T)
        outx = STC(mu = mu, newcombinations = out$combinations, priority = priority,disjointcom = disjointcom,
                   combinations = combinations, group = NULL, plus = F)
        name = out$name
        name.joint = out$name.joint
        joint = out$joint
        newcombinations = out$combinations
        newcombinationsplus = out1$newcombinations
        newcombinationsminus = outx$newcombinations
        if(all(newcombinationsplus == newcombinationsminus)){eq = 0}else{eq = 1}
        newmc = out1$newmc
        mc = out1$mc
        combinationsgroup = out1$combinationsgroup
        #combinations
        combinationsplus = out1$combinations
        combinationsminus = outx$combinations
      }
    }else{newcombinations = combinations;newmc = mc}
    for(i in 1:length(joint)){
      name1 = str_split(name.joint[i],"")[[1]]
      value_a = disjointcom[which(names(disjointcom)==name1[1])]
      value_b = disjointcom[which(names(disjointcom)==name1[2])]
      if(value_a<joint[i] || value_b<joint[i]){stop("joint part cannot be larger than the whole part")}
    }
    group = RCH(name.joint = name.joint, name.dis= name.dis)
    grouplength = do.call(c,lapply(group, length))
    group = group[order(grouplength,decreasing = T)]
    grouplength = do.call(c,lapply(group, length))
    groupnum = length(group)
    groupxy = list()
    radiusxy = list()
    namexy = list()
    area = areacompute(newcombinations = newcombinations, disjointcom = disjointcom, ThreeD = ThreeD,
                       newmc = newmc, fit = fit)
    if(length(priority) == 1 && priority==2){
      combinationsgroup = lapply(group,function(a,mu,name.dis,name){y = name.dis[a];for(i in 1:length(name)){
        if(any(y%in%str_split(name[i],"")[[1]])){y = c(y,name[i])}};y = unique(y);z = str_length(y);
        y = combinations[which(name%in%y)];return(y)},mu = mu,name.dis = name.dis, name = names(combinations))
      NH = list()
      fake = rep(0,40)
      for(fitnumber in 1:40){
        out2 = DC(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,
                  newcombinations = newcombinations,NAME = NAME, ALPHA = ALPHA,
                  radius = area$radius, ThreeD = ThreeD, joint.area = area$joint.area, Area = area$Area,
                  GRAM = GRAM, arbitrary = arbitrary, fitnumber = fitnumber)
        lossnum = out2$lossnum
        fake[fitnumber] = lossnum
        NH[[fitnumber]] = out2
        if(fitnumber == 1){next}
        else if(fake[fitnumber] >= fake[fitnumber-1]){out2 = NH[[fitnumber-1]];break}
        else{out2 = NH[[fitnumber]]}
      }
      if(priority%in%mc){
        weighted.least.square = HH(combinations = combinations, disjointcom = disjointcom, priority = priority,
                                   combinationsgroup = combinationsgroup, out2 = out2, weight = weight)
      }else{weighted.least.square = 0}
    }else{
      out4 = PAS(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,priority = priority,
                 weight = weight, newcombinations = newcombinations,NAME = NAME, ALPHA = ALPHA,fit = fit,area = area,ThreeD = ThreeD,
                 GRAM = GRAM, arbitrary = arbitrary, newcombinationsplus = newcombinationsplus, mu = mu, eq = eq,combinations = combinations,
                 newcombinationsminus = newcombinationsminus,disjointcom = disjointcom, newmc = newmc,
                 combinationsgroup = combinationsgroup, combinationsplus = combinationsplus, combinationsminus = combinationsminus)
      out2 = out4$out2
      weighted.least.square = out4$weighted.least.square
    }
    lossnum = out2$lossnum
    groupxy = out2$groupxy
    radiusxy = out2$radiusxy
    namexy = out2$namexy
    #combine all groups
    groupxylist = groupxy
    radiusxylist = radiusxy
    if(groupnum>1){
      out3 = QNC(groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy,radius = area$radius, delta = delta, namexy = namexy)
      groupxylist = out3$groupxylist
      radiusxylist = out3$radiusxylist
      namexy = out3$namexy
    }
    xy = groupxylist[[groupnum]]
    radius = radiusxylist[[groupnum]]
    name.dis = namexy[[groupnum]]
    if(ThreeD){
      open3d()
      spheres3d(xy[,1],xy[,2],xy[,3], radius = radius,
                color = cols, alpha = alpha.col)
      text3d(xy[,1],xy[,2],xy[,3],name.dis)
    }else{
      plot.new()
      circle.plot(xy,radius,name.dis,cols, main = main, mar = mar)
    }
    outfinal = list(center = xy,diameter = radius, LOSS = lossnum, weighted.least.square = unname(weighted.least.square))
  }
  return(outfinal)
}

#--- helper functions -----------------------------------------------------

#library(rgl)
#library(strinr)
#library(Rcpp)
#sourceCpp(file = "your path/vennCpp2.cpp")

#Big loop
DC = function(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,
              newcombinations = newcombinations, NAME = NAME, ALPHA = ALPHA,radius = radius, ThreeD = ThreeD,
              joint.area = joint.area, Area = Area,GRAM = GRAM, arbitrary = arbitrary, fitnumber = fitnumber){
  mc = str_length(names(newcombinations))
  name.dis = names(newcombinations)[which(mc==1)]
  name.joint = names(newcombinations)[which(mc!=1)]
  lossnum = 0
  for(gn in 1:groupnum){
    name.dis1 = name.dis[group[[gn]]]
    m1 = length(name.dis1)
    if(m1 == 1){
      if (is.null(NAME)){name1 = name.dis1}else{name1 = NAME[group[[gn]]]}
      radius1 = radius[name.dis1]
      if(ThreeD){groupxy[[gn]] = matrix(c(0,0,0),nrow = 1)}else{groupxy[[gn]] = matrix(c(0,0),nrow=1)}
      radiusxy[[gn]] = radius1
      namexy[[gn]] = name1
      rownames(groupxy[[gn]]) = name1
    }else{
      name.joint1 = rep(0,length(name.joint))
      for (i in 1:length(name.joint)){
        if(any(str_detect(name.joint[i],name.dis1))){
          name.joint1[i] = name.joint[i]
        }
      }
      joint.area1 = joint.area[which(names(joint.area) == name.joint1)]
      if(length(which(name.joint1 == "0"))!=0){
        name.joint1 = name.joint1[-which(name.joint1 == "0")]
      }
      Area1 = Area[name.dis1]
      radius1 = radius[name.dis1]
      N = matrix(rep(0,length(name.joint1)*length(name.dis1)),nrow = length(name.joint1))
      for (i in 1:length(name.joint1)){
        for (j in 1:length(name.dis1)){
          if(str_detect(name.joint1[i],name.dis1[j])){
            N[i,j] = j
          }
        }
      }
      N=t(N)
      N = matrix(N[which(N != 0)],ncol = 2, byrow = T)
      ED = matrix(rep(0, m1^2),ncol = m1)
      W = matrix(rep(0,m1^2),nrow = m1,ncol = m1)
      D = W
      for (k in 1:dim(N)[1]){
        for (i in 1:m1){
          for (j in 1:m1){
            if (i == N[k,1] && j == N[k,2]){
              #Euler distance
              ED[i,j] = Distance(radius1[i],radius1[j], joint.area1[k],ThreeD = ThreeD)
              #Jaccard distance
              D[i,j] = joint.area1[k]/(Area1[i]+Area1[j]-joint.area1[k])
            }
          }
        }
      }
      for (i in 1:(m1-1)){
        for (j in (i+1):m1){
          if (ED[i,j] == 0){
            ED[i,j] = Distance(radius1[i],radius1[j],0,ThreeD = ThreeD)
          }
          else{next}
        }
      }
      ED = ED+t(ED)
      ED[which(ED==10)]=0
      rownames(ED) = name.dis1
      colnames(ED) = name.dis1
      D = D+t(D)+diag(1,m1)
      rownames(D) = name.dis1
      colnames(D) = name.dis1
      D = 1-D
      #Using Gram matrix to find the initial location
      if (GRAM){
        D2 = ED^2
        J = diag(m1)-1/m1
        G = -0.5*J%*%D2%*%J
        Em = svd(G)$u
        lambdam = svd(G)$d
        U = Em*sqrt(lambdam)
      }else{
        if(is.null(arbitrary)){
          set.seed("12345")
          kk = floor(runif(1,1,m1+1))
        }else{
          kk = arbitrary
        }
        for (i in 1:m1){
          for (j in 1:m1){
            W[i,j] = (D[i,kk]^2+D[j,kk]^2-D[i,j]^2)/2
          }
        }
        U = svd(W)$u
      }
      #Standardize coordinates
      if(ThreeD){
        if (dim(U)[1]>=3){
          xy = U[,1:3]
        }else{
          xy = cbind(U,rep(0,2))
        }
      }else{
        xy = U[,1:2]
        x = scale(xy[,1])
        y = scale(xy[,2])
        a1 = min(min(x-fitnumber*radius1),min(y-fitnumber*radius1))
        a2 = max(max(x+fitnumber*radius1),max(y+fitnumber*radius1))
        #Standardize Scale in (0,1) vs (0,1)
        x = (x-a1)/(a2-a1)
        y = (y-a1)/(a2-a1)
        xy = cbind(x,y)
      }
      L = LoopR(xy = xy, ALPHA = ALPHA, radius = radius1, ED = ED, ThreeD = ThreeD)
      if (is.null(NAME)){name1 = name.dis1}else{name1 = NAME[group[[gn]]]}
      groupxy[[gn]] = L$xy
      radiusxy[[gn]] = radius1
      namexy[[gn]] = name1
      rownames(groupxy[[gn]]) = name1
      lossnum = lossnum + L$f1
    }
  }
  out = list(lossnum = lossnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy)
}

MC = function(combinations = combinations,disjointcom = disjointcom, mu = mu){
  mc = str_length(names(combinations))
  joint2 = combinations[which(mc == 2)]
  name.2 = names(joint2)
  m3 = unique(mc)
  m3 = m3[which(m3>2)]
  resam = 0
  for (i in 1:length(m3)){
    new_joint = matrix(rep(0,10000),nrow = 2)
    jointn = sort(combinations[which(mc == m3[i])],decreasing = T)
    com = combn(m3[i],2)
    L = dim(com)[2]
    for (j in 1:length(jointn)){
      for (k in 1:L){
        name1 = str_split(names(jointn[j]),"")[[1]][com[,k]]
        value_a = disjointcom[which(names(disjointcom)==name1[1])]
        value_b = disjointcom[which(names(disjointcom)==name1[2])]
        valuemin = min(value_a,value_b)
        name2 = paste(name1,collapse="")
        if(length(which(name.2==name2))==0){
          resam = resam+1
          if(jointn[j] == 0){
            jointn[j] = min(disjointcom)*0.01
          }
          if((jointn[j]*mu^(m3[i]-2)) < valuemin){
            new_joint[1,resam] = (jointn[j]*mu^(m3[i]-2))
          }
          else if(jointn[j]>valuemin){
            stop("joint part cannot be larger than the whole part")
          }
          else{
            new_joint[1,resam] = runif(1,jointn[j],valuemin)
          }
          new_joint[2,resam] = name2
        }
        else if(joint2[which(name.2==name2)]<jointn[j]){
          stop("higher way intersection cannot be larger than lower way one")
        }
      }
    }

    new_joint = new_joint[which(new_joint!=0)]
    joint_new = rep(0,length(new_joint)/2)
    name.joint_new = rep(0,length(new_joint)/2)
    for (i in 1:length(joint_new)){
      joint_new[i] = as.numeric(new_joint[2*i-1])
      name.joint_new[i] = new_joint[2*i]
    }
    names(joint_new) = name.joint_new
    joint2 = c(joint2,joint_new)
    name.2 = c(name.2,name.joint_new)
  }
  combinations = c(disjointcom,joint2)
  name = names(combinations)
  mc = str_length(name)
  name.joint = name[-which(mc==1)]
  joint = combinations[-which(mc==1)]
  ######################
  if(any(is.na(name))){name = name[-which(is.na(name))]}
  if(any(is.na(name.joint))){name.joint = name.joint[-which(is.na(name.joint))]}
  if(any(is.na(joint))){joint = joint[-which(is.na(joint))]}
  if(any(is.na(combinations))){combinations = combinations[-which(is.na(combinations))]}
  if(any(is.na(mc))){mc = mc[-which(is.na(mc))]}

  out = list(name= name, name.joint = name.joint, joint = joint, combinations = combinations,mc = mc,resam = resam)
  return(out)
}


#transform disjoint
disjoint.transform = function(combinations = NULL){
  name = names(combinations)
  mc = str_length(name)
  combinations = combinations[order(mc)]
  name = names(combinations)
  mc = sort(mc)
  if(max(mc)!=1){
    for(i in 1:length(mc)){
      if(mc[i] == max(mc)){break}else{
        ma = sapply(str_split(names(combinations[-i]),pattern = ""),
                    function(a){y = 0; if(all(str_split(names(combinations[i]),pattern = "")[[1]]%in%a)){y = 1};return(y)})
        if(length(which(ma == 1))!=0){combinations[i] = combinations[i] + sum(combinations[-i][which(ma == 1)])}
      }
    }
  }
  return(combinations)
}

# If input is a data list
datatocom = function(data = NULL){
  com = t(combn(dim(data)[2],2))
  letter = c(LETTERS,letters)
  M = rep(0, dim(com)[1])
  Names = rep(0, dim(com)[1])
  for(j in 1:dim(com)[1]){
    n = 0
    for(i in 1:dim(data)[1]){
      if(data[i,com[j,][1]] == 1 && data[i,com[j,][2]] == 1){
        n = n+1
      }
    }
    M[j] = n
    Names[j] = paste(c(letter[com[j,][1]],letter[com[j,][2]]),collapse="")
  }
  names(M) = Names
  tot = apply(data,2,"sum")
  names(tot) = letter[1:dim(data)[2]]
  combinations = c(tot, M)
  if(is.null(names(data))){
    Name = Names
  }else{
    Name = names(data)
  }

  out = list(combinations = combinations, Name = Name)
  return(out)
}

# Calculates distance between two circles
Distance = function(r1,r2,S,ThreeD = ThreeD){
  theta1 = seq(0,pi, length = 100)
  theta2 = seq(0,pi, length = 100)
  if (ThreeD){
    f1 = matrix(rep(pi/3*r2^3*(1-cos(theta2))^2*(2+cos(theta2)),100), nrow = 100, byrow = T) +
      pi/3*r1^3*(1-cos(theta1))^2*(2+cos(theta1)) - S
    f2 = matrix(rep(r1*sin(theta1),100), nrow = 100, byrow = T) - r2*sin(theta2)
  }
  else{
    f1 = matrix(rep(theta2*r2^2 - sin(2*theta2)*r2^2/2,100), nrow = 100, byrow = T) +
      theta1*r1^2 - sin(2*theta1)*r1^2/2 - S
    f2 = matrix(rep(r1*sin(theta1),100), nrow = 100, byrow = T) - r2*sin(theta2)
  }
  R = abs(f1) + abs(f2)
  k = which(R== min(R),arr.ind=T)[1,]
  d = r1*cos( theta1[k[1]]) + r2*cos(theta2[k[2]])
  if (d==0){d=10}
  return(d)
}

# Separate each group
M3 = function(i = i,clustermatrix = clustermatrix){
  index = which(clustermatrix[i,]!= 0)
  indexi = rep(0,length(index)+1)
  while(length(index) != length(indexi)){
    for(j in 1:length(index)){
      if(j==1){indexi = index}
      indexi = unique(c(indexi, which(clustermatrix[index[j],]!= 0)))
    }
    if(length(index) == length(indexi)){break}
    else{
      index=indexi
      indexi = rep(0,length(index)+1)}
  }
  return(index)
}
RCH = function(name.joint, name.dis){
  m = length(name.dis)
  clustermatrix = matrix(rep(0,m*m),nrow = m)
  for(i in 1:m){
    namei = name.dis[i]
    loca = unique(str_split(paste(name.joint[which(str_detect(name.joint,namei)==T)],
                                  collapse=""),"")[[1]])
    loca1 = rep(0,length(loca))
    for(j in 1:length(loca)){
      loca1[j] = which(name.dis == loca[j])
    }
    clustermatrix[i,loca1] = 1
    if(clustermatrix[i,i] == 0){
      clustermatrix[i,i] = 1
    }
  }
  group = list()
  for(i in 1:m){
    if(i ==1){
      index = M3(i, clustermatrix)
      group[[i]] = sort(index)
      ID = index
      groupnum = i
    }
    else if(length(which(ID == i)) == 0){
      index = M3(i, clustermatrix)
      group[[groupnum+1]] = sort(index)
      ID = c(ID,index)
      groupnum = groupnum+1
    }
  }
  return(group)
}

# Plots circles of Venn diagram
circle.plot =  function(xy = NULL, radius =NULL,name = NULL, col = NULL, main=NULL, mar = NULL) {
  par(mar = mar)
  a1 =  range(c((xy + radius)[,1], (xy - radius)[,1]))
  a2 = range(c((xy + radius)[,2], (xy - radius)[,2]))
  plot.window(a1, a2, "", asp = 1)
  theta = seq.int(360)/360*2*pi
  for (i in 1:length(radius)){
    polygon(xy[i,1] +  radius[i]*cos(theta), xy[i,2] + radius[i]*sin(theta), col = col[i],border = col[i])
  }
  text(xy, name)
  title(main = main)
}

Rotate = function(xy,transxy,radius1,radius2, delta){
  for(i in 1:35){
    theta = i/18*pi
    if(dim(xy)[2] == 3){
      rotation1 = matrix(c(1,0,0,0,cos(theta),-sin(theta),0,sin(theta),cos(theta)),nrow = 3,byrow = T)
      rotation2 = matrix(c(cos(theta),0,sin(theta),0,1,0,-sin(theta),0,cos(theta)),nrow = 3,byrow = T)
      rotation3 = matrix(c(cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1),nrow = 3,byrow = T)
      rotation = rotation1%*%rotation2%*%rotation3}
    else{
      rotation = matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow = 2,byrow = T)
    }
    newxy = t(rotation%*%t(transxy))
    newxy = t(transxy[1,] - newxy[1,] + t(newxy))
    out = alldisR(xy,newxy,radius1,radius2, delta)
    if (out!=0 ){break}
  }
  la = list(out =out, newxy = newxy, i =i)
  return(la)
}

EV3 = function(combinations = combinations, m =m, ThreeD = ThreeD, cols = cols, main = main, mar = mar,
               fit = fit, name.dis = name.dis,delta = delta, NAME=NAME, alpha.col = alpha.col){
  #Compute a and c
  mc = str_length(names(combinations))
  disjointcom = combinations[which(mc==1)]
  a = combinations
  for (j in 1:length(combinations)){
    a[j] = combinations[j]/(sum(disjointcom))
  }
  #Get each circle's radius
  if (ThreeD){
    Area = a[which(mc==1)]
    radius = (3*Area/(4*pi))^(1/3)
  }else{
    Area = a[which(mc==1)]
    radius = sqrt(Area/pi)*fit
    a = a*fit^2
  }
  if(ThreeD){
    xy = matrix(rep(0,3),nrow = 1)
    xy1 = matrix(rep(0,3),nrow = 1)
  }else{
    xy = matrix(rep(0,2),nrow = 1)
    xy1 = matrix(rep(0,2),nrow = 1)
  }
  for(i in 1:(m-1)){
    radiusvec = radius[i+1]
    xyvec = transR(xy = xy, radius = radius[1:i], radiusvec = radiusvec, radiusall = radius)
    xy = rbind(xy,xyvec)
  }
  rownames(xy) = name.dis
  for(i in 1:(m-1)){
    radiusvec = radius[i+1]
    direc = -xy[i+1,]
    direc = direc/sqrt(sum(direc^2))*delta/10
    newxy = closeR(xy1,matrix(xy[i+1,],nrow=1),radius[1:i],radius[i+1],delta = delta,direc)$xy
    xy1 = rbind(xy1,newxy)
  }
  if (is.null(NAME)){name = name.dis}else{name = NAME}
  if(ThreeD){
    open3d()
    spheres3d(xy1[,1],xy1[,2],xy1[,3], radius = radius,
              color = cols, alpha = alpha.col)
    text3d(xy1[,1],xy1[,2],xy1[,3], name)
  }
  else{
    plot.new()
    circle.plot(xy1,radius,name = name, col = cols, main = main, mar = mar)
  }
  out = list(xy = xy1, radius = radius)
  return(out)
}


QNC = function(groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy,radius = radius, delta = delta, namexy = namexy){
  groupxylist = groupxy
  radiusxylist = radiusxy
  for (gn in 1:(groupnum-1)){
    minx = which(groupxy[[gn+1]][,1] == min(groupxy[[gn+1]][,1]))
    radiusvec = radiusxy[[gn+1]][minx]
    xyvec = transR(xy = groupxy[[gn]], radius = radiusxy[[gn]], radiusvec = radiusvec, radiusall = radius)
    #move
    transxy = t(t(groupxy[[gn+1]]) + xyvec - groupxy[[gn+1]][minx,])
    if(alldisR(groupxy[[gn]],transxy,radiusxy[[gn]],radiusxy[[gn+1]],delta) == 0){
      out = Rotate(groupxy[[gn]],transxy,radius[[gn]],radius[[gn+1]],delta)
      outr = out$out
      newxy = out$newxy
      while(outr == 0){
        xyvec = transR(xy = groupxy[[gn]], radius = radiusxy[[gn]], radiusvec = radiusvec, radiusall = radius)
        #move
        transxy = t(t(groupxy[[gn+1]]) + xyvec - groupxy[[gn+1]][minx,])
        out = Rotate(groupxy[[gn]],transxy,radius[[gn]],radius[[gn+1]],delta)
        outr = out$out
        newxy = out$newxy
      }
      rownames(newxy) = namexy[[gn+1]]
    }else{newxy = transxy}
    groupxy[[gn+1]] = rbind(groupxy[[gn]],newxy)
    groupxylist[[gn+1]] = newxy
    radiusxy[[gn+1]] = c(radiusxy[[gn]],radiusxy[[gn+1]])
  }
  for(gn in 1:(groupnum-1)){
    center1 = rep(0,dim(groupxylist[[gn]])[2])
    center2 = rep(0,dim(groupxylist[[gn+1]])[2])
    for(j in 1:length(center1)){
      center1[j] = (max(groupxylist[[gn]][,j]+radiusxylist[[gn]])+min(groupxylist[[gn]][,j]-radiusxylist[[gn]]))/2
      center2[j] = (max(groupxylist[[gn+1]][,j]+radiusxylist[[gn+1]])+min(groupxylist[[gn+1]][,j]-radiusxylist[[gn+1]]))/2
    }
    direc = center1-center2
    if(delta==0){delta=10^(-4)}
    direc = direc/sqrt(sum(direc^2))*delta/10
    newxy = closeR(groupxylist[[gn]],groupxylist[[gn+1]],radiusxylist[[gn]],radiusxylist[[gn+1]],delta = delta,direc)$xy
    groupxylist[[gn+1]] = rbind(groupxylist[[gn]],newxy)
    radiusxylist[[gn+1]] = c(radiusxylist[[gn]],radiusxylist[[gn+1]])
    namexy[[gn+1]] = c(namexy[[gn]],namexy[[gn+1]])
  }
  out = list(groupxylist = groupxylist, radiusxylist = radiusxylist,namexy = namexy)
  return(out)
}

DP = function(y = y, z=z, mu = mu, plus = plus){
  name2 = names(y[which(z==2)])
  n3 = y[which(z>2)]
  if(length(name2) != 0 && length(n3) != 0){
    n1 = y[which(z==1)]
    n2 = NULL
    name3 = names(n3)
    for(i in 1:length(name3)){
      x = rep(0,length(name2))
      for(j in 1:length(name2)){
        if(all(str_split(name2[j],pattern = "")[[1]]%in%str_split(name3[i],pattern = "")[[1]])){
          x[j] = j
        }
      }
      x = x[which(x!=0)]
      if(length(x) == dim(combn(str_length(name3[i]),2))[2] && plus == TRUE){
        nx = y[which(z==2)][x] + (mu-0.1)*n3[i]
        if(any(nx>min(n1[which(names(n1)%in%str_split(names(nx),pattern = "")[[1]]==T)]))){
          nx = y[which(z==2)][x]
        }
        n2 = c(n2,nx)
      }else if(length(x) == dim(combn(str_length(name3[i]),2))[2] && plus == FALSE){
        nx = y[which(z==2)][x] - (mu-0.1)*n3[i]
        if(any(nx<n3[i])){
          nx = y[which(z==2)][x]
        }
        n2 = c(n2,nx)
      }
    }
    n2 = c(n2,y[which(z==2)][which(names(y[which(z==2)])%in%names(n2)==F)])
    if(is.null(n2) == FALSE){
      y = c(n1,n2,n3)
    }
  }
  return(y)
}

#combinations separate
STC = function(mu = mu, newcombinations = newcombinations, priority = priority,disjointcom = disjointcom,
               combinations = combinations, group = NULL, plus = TRUE){
  sepgroup = function(a,x=0, mu = mu, plus = plus){
    y = name.dis[a]
    for(i in 1:length(name)){
      if(any(y%in%str_split(name[i],"")[[1]])){
        y = c(y,name[i])
      }
    }
    y = unique(y)
    z = str_length(y)
    y = combinations[which(name%in%y)]
    if(x==0){y = DP(y=y,z=z, mu = mu, plus = plus)}
    return(y)
  }
  mc = str_length(names(newcombinations))
  name.dis = names(newcombinations)[which(mc==1)]
  name.joint = names(newcombinations)[which(mc!=1)]
  name = names(combinations)
  if(is.null(group)){
    group = RCH(name.joint = name.joint, name.dis= name.dis)
    grouplength = do.call(c,lapply(group, length))
    group = group[order(grouplength,decreasing = T)]
    grouplength = do.call(c,lapply(group, length))
    groupnum = length(group)
    combinationsgroup = lapply(group, sepgroup, x = 1, mu = mu, plus = plus)
    combinations = unlist(lapply(group, sepgroup,x = 0, mu = mu, plus = plus))
    mc = str_length(names(combinations))
    combinations = combinations[order(mc)]
    mc = sort(mc)
    out = MC(combinations = combinations,disjointcom = disjointcom,mu = mu)
    name = out$name
    name.joint = out$name.joint
    joint = out$joint
    newcombinations = out$combinations
    newmc = out$mc
    co = list(name = name, name.joint = name.joint, joint = joint, newcombinations = newcombinations, newmc = newmc,
              combinationsgroup = combinationsgroup, combinations = combinations, mc = mc)
  }else{
    combinations = unlist(lapply(group, sepgroup,x = 0, mu = mu, plus = plus))
    mc = str_length(names(combinations))
    combinations = combinations[order(mc)]
    mc = sort(mc)
    out = MC(combinations = combinations,disjointcom = disjointcom,mu = mu)
    name = out$name
    name.joint = out$name.joint
    joint = out$joint
    newcombinations = out$combinations
    newmc = out$mc
    co = list(name = name, name.joint = name.joint, newcombinations = newcombinations, newmc = newmc,
              combinations = combinations, mc = mc)
  }
  return(co)
}

pixel = function(xy = xy,radius = radius,num=100,combinations = combinations, weight = weight,priority = priority,
                 name.priorityway = name.priorityway, i = i, priorityway = priorityway){
  name = names(combinations)
  name.dis = name[which(str_length(names(combinations))==1)]
  m = length(name.dis)
  l = list()
  xuan = (max(xy[,1]) - min(xy[,1])+2*max(radius))/num
  yuan = (max(xy[,2]) - min(xy[,2])+2*max(radius))/num
  for (k in 1:m){
    l[[k]] = listR1(matrix(rep(0,num^2),nrow =num),xy,radius,k,yuan,xuan,num=num)
  }
  M = listR2(myList =  l, m = m, num = num)
  M = M[which(listR4(M)!=0),]
  Me = matrix(unlist(unique(as.data.frame(M))),ncol = m)
  Len_Me = listR3(M,Me)
  names(Len_Me) = apply(Me,1,function(a){paste(name.dis[which(a!=0)],collapse="")})
  Len_Me = Len_Me/sum(Len_Me)
  mk = Len_Me[which(str_length(names(Len_Me))>=min(str_length(name.priorityway[[i]])))]
  lo = 0
  if(length(mk) !=0){
    for(j in 1:length(mk)){
      #weight
      mb = which(priority == str_length(names(mk[j])))
      if(length(mb)==1){
        w = weight[mb]
      }else{w = 1}
      #loss
      ma = sapply(str_split(names(mk[-j]),pattern = ""),
                  function(a){y = 0; if(all(str_split(names(mk[j]),pattern = "")[[1]]%in%a)){y = 1};return(y)})
      if(length(which(ma == 1))!=0){mk[j] = mk[j] + sum(mk[-j][which(ma == 1)])}
      mb = which(names(priorityway[[i]])==names(mk[j]))
      if(length(mb)!=0){
        lo = lo + (priorityway[[i]][mb]-mk[j])^2*w
      }else{
        lo = lo + mk[j]^2*w
      }
    }
  }else{
    for(j in 1:length(priorityway[[i]])){
      lo = lo + priorityway[[i]][j]^2*weight[which(priority%in%str_length(names(priorityway[[i]][j])))]
    }
  }
  return(lo)
}

HH = function(combinations = combinations, disjointcom = disjointcom, priority = priority,
              combinationsgroup = combinationsgroup, out2 = out2, weight = weight){
  mc = str_length(names(combinations))
  Losspriority = 0
  truea = lapply(combinationsgroup,function(a){a/sum(disjointcom)})
  priorityway = lapply(truea, function(a){a[which(str_length(names(a))%in%priority==TRUE)]})
  name.priorityway =  lapply(priorityway, function(a){names(a)})
  groupin = sapply(priorityway, function(a){y = 0;if(length(a)!=0){y = 1};return(y)})
  nonzero = which(groupin!=0)
  groupinxy = list()
  for(i in 1:length(nonzero)){
    Losspriority = Losspriority + pixel(xy = out2$groupxy[[nonzero[i]]], radius = out2$radiusxy[[nonzero[i]]], num = 100,
                                        combinations = combinationsgroup[[nonzero[i]]], weight = weight,priority = priority,
                                        name.priorityway = name.priorityway, i = i, priorityway = priorityway)
  }
  return(Losspriority)
}

areacompute = function(newcombinations = newcombinations, disjointcom = disjointcom,
                       ThreeD = ThreeD, newmc = newmc, fit = fit){
  a = newcombinations/sum(disjointcom)
  #Get each circle's radius
  if (ThreeD){
    Area = a[which(newmc==1)]
    radius = (3*Area/(4*pi))^(1/3)
  }else{
    Area = a[which(newmc==1)]
    radius = sqrt(Area/pi)*fit
    a = a*fit^2
  }
  joint.area = a[which(newmc!=1)]
  Area = a[which(newmc==1)]
  out = list(radius = radius, a = a, joint.area = joint.area, Area = Area)
  return(out)
}


PAS = function(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,priority = priority,
               weight = weight,newcombinations = newcombinations,NAME = NAME, ALPHA = ALPHA,fit = fit,area = area,ThreeD = ThreeD,
               GRAM = GRAM, arbitrary = arbitrary, newcombinationsplus = newcombinationsplus, mu = mu, eq = eq,combinations = combinations,
               newcombinationsminus = newcombinationsminus,disjointcom = disjointcom, newmc = newmc,
               combinationsgroup = combinationsgroup, combinationsplus = combinationsplus, combinationsminus = combinationsminus){
  NH = list()
  fake1 = rep(0,40)
  for(tu in 1:40){
    if(tu ==1){
      NH1 = list()
      fake = rep(0,40)
      for(fitnumber in 1:40){
        out2 = DC(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,
                  newcombinations = newcombinations,NAME = NAME, ALPHA = ALPHA,
                  radius = area$radius, ThreeD = ThreeD, joint.area = area$joint.area, Area = area$Area,
                  GRAM = GRAM, arbitrary = arbitrary, fitnumber = fitnumber)
        lossnum = out2$lossnum
        fake[fitnumber] = lossnum
        NH1[[fitnumber]] = out2
        if(fitnumber == 1){next}
        else if(fake[fitnumber] >= fake[fitnumber-1]){out2 = NH1[[fitnumber-1]];break}
        else{out2 = NH1[[fitnumber]];fitnumber = fitnumber+1}
      }
      #search deriction
      if(eq == 1){
        area = areacompute(newcombinations = newcombinationsplus, disjointcom = disjointcom, ThreeD = ThreeD,
                           newmc = newmc, fit = fit)
        out2plus = DC(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,
                      newcombinations = newcombinationsplus, NAME = NAME, ALPHA = ALPHA,
                      radius = area$radius, ThreeD = ThreeD, joint.area = area$joint.area, Area = area$Area,
                      GRAM = GRAM, arbitrary = arbitrary, fitnumber = fitnumber-1)
        area = areacompute(newcombinations = newcombinationsminus, disjointcom = disjointcom, ThreeD = ThreeD,
                           newmc = newmc, fit = fit)
        out2minus = DC(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,
                       newcombinations = newcombinationsminus, NAME = NAME, ALPHA = ALPHA,
                       radius = area$radius, ThreeD = ThreeD, joint.area = area$joint.area, Area = area$Area,
                       GRAM = GRAM, arbitrary = arbitrary, fitnumber = fitnumber-1)
        ff1 = HH(combinations = combinations, disjointcom = disjointcom, priority = priority,
                 combinationsgroup = combinationsgroup, out2 = out2plus, weight = weight)
        ff2 = HH(combinations = combinations, disjointcom = disjointcom, priority = priority,
                 combinationsgroup = combinationsgroup, out2 = out2minus, weight = weight)
        ff3 = HH(combinations = combinations, disjointcom = disjointcom, priority = priority,
                 combinationsgroup = combinationsgroup, out2 = out2, weight = weight)
        if(ff3<=ff1 && ff3<=ff2){fake1[tu] = ff3;break}else if(ff1<ff2){combinations1 = combinationsplus;plus = T;
        newcombinations = newcombinationsplus;fake1[tu] = ff1;
        }else{combinations1 = combinationsminus;plus = F;newcombinations = newcombinationsminus;fake1[tu] = ff2}
      }else{fake1[tu] = HH(combinations = combinations, disjointcom = disjointcom, priority = priority,
                           combinationsgroup = combinationsgroup, out2 = out2, weight = weight)};
      NH[[tu]] = out2
    }else{
      mu = mu + 0.1
      if(eq == 1){
        out1 = STC(mu = mu, newcombinations = newcombinations, priority = priority,disjointcom = disjointcom,
                   combinations = combinations1, group = group, plus = plus)
        combinations1 = out1$combinations
        newcombinations = out1$newcombinations
        newmc = out1$newmc
      }else{
        out = MC(combinations = combinations,disjointcom = disjointcom,mu = mu)
        newcombinations = out$combinations
        newmc = out$mc
      }
      area = areacompute(newcombinations = newcombinations, disjointcom = disjointcom,
                         ThreeD = ThreeD, newmc = newmc, fit = fit)
      out2 = DC(group = group, groupnum = groupnum, groupxy = groupxy, radiusxy = radiusxy, namexy = namexy,
                newcombinations = newcombinations,NAME = NAME, ALPHA = ALPHA,
                radius = area$radius, ThreeD = ThreeD, joint.area = area$joint.area, Area = area$Area,
                GRAM = GRAM, arbitrary = arbitrary, fitnumber = fitnumber - 1)
      fake1[tu] = HH(combinations = combinations, disjointcom = disjointcom, priority = priority,
                     combinationsgroup = combinationsgroup, out2 = out2, weight = weight)
    }
    NH[[tu]] = out2
    if(tu==1 || tu==2){next}
    else if(fake1[tu] >= fake1[tu-1]){out2 = NH[[tu-1]];break}else{out2 = NH[[tu]]}
  }
  if(tu == 1){weighted.least.square = fake1[tu]}else{weighted.least.square = fake1[tu-1]}
  out = list(weighted.least.square = weighted.least.square, out2 = out2)
  return(out)
}
