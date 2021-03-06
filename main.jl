# Rounded Edge Airfoil generator
# Developed by Yi Tsung Lee, NCSU AE
# The generated file will have the LE coordinate near (0,0)
# The generated file will have the TE coordinate at (c,0)
# outputfile will be in xfoil format(out_type =1)
# outputfile include num point in the section and 3-d coordinate (out_type =2)

# Remenber to check the working directory first for file saving

#=  In case these Pkg not existed, activate this first
Pkg.add("Plots")
Pkg.add("QuadGK")   # Numerical intergration Pkg
Pkg.add("Calculus") # Derrivative Pkg
=#

using Plots
using DelimitedFiles
using LinearAlgebra
#using QuadGK
#using Calculus
gr()

function main()

   # =============== control variable section ================
   # out_type = 1:xfoil data form; = 2:3_section data
   out_type = 1

   # variable definition
   c::Float64      = 1.0 # chord length
   max_t::Float64  = 0.12 # max thickness (between 0~1.0)
   max_tp::Float64 = 0.0 # max thickness position(in x/c,should be 0~1.0)
   num_p::Int64    = 200 # number of points
   camber::Float64 = 0.06 # camber (in y/c )

   # =============== control variable section ================
   # start of the coordinate generation

   ellipse(num_p,c,max_t,camber,out_type)

end

# elliptical shape generator
function ellipse(num_p::Int64,c::Float64,max_t::Float64,
                 camber::Float64,out_type::Int64)
 # initialize
  num_up = round(Int64,num_p/2+1)
  num_lp = num_p-num_up+2

  theta = zeros(Float64,num_p+1)  # theta used to generate the x point value
  xy  = zeros(Float64,num_p+1,2)  # xy data set (for xfoil),num_p+1 is identical to point 1
  xy_c = zeros(Float64,num_p+1,2) # xy data for camber line
  xyz_l = zeros(Float64,num_lp,3) # xyz data set for lower surface
  xyz_u = zeros(Float64,num_up,3) # xyz data set for upper surface

 # ellipse based coordinate
  # equal spacing theta generation
 x_ellip(a) = c/2*(1+cos(a)) # elliptical x function (a= theta here)
 y_ellip(a) = max_t/2*sin(a) # elliptical y function (a= theta here)
 r_p = 0.25  # refinement portion, should be between 0~25%
 r_ratio = refinement_iteration(x_ellip,y_ellip,num_p,r_p)

  # initialize the theta
  for i  = 1:num_p+1
     theta[i] = 2*pi/num_p*(i-1)
  end
  #  do the actual refinement
  theta = refinement(r_ratio,1,round(Int64,num_p*r_p),theta)
  theta = refinement(r_ratio,round(Int64,num_p*0.5),round(Int64,num_p*(0.5-r_p)),theta)
  theta = refinement(r_ratio,round(Int64,num_p*0.5),round(Int64,num_p*(0.5+r_p)),theta)
  theta = refinement(r_ratio,num_p+1,round(Int64,num_p*(1-r_p)),theta)

  for i  = 1:num_p+1
      xy[i,1] = x_ellip(theta[i])
      xy[i,2] = y_ellip(theta[i])
  end

  # apply parabolic camberline

  # analytical intergration
  k = 4*camber/c
  y_para(x) = -k*x*(x-c)  # parabolic equation function
  dy_para(x) = -2*k*x+c*k   # parabolic slope-2*k+c*k function

  # store the expression of s(x) and its differentiate function
  s(x) = 1/2*x*sqrt(1+(2*k*x)^2)+1/(4*k)*log(sqrt(1+(2*k*x)^2)+2*k*x)
  ds(x) = 1/2*sqrt(1+(2*k*x)^2) + 1/2*x*(4*k^2*x/sqrt(1+(2*k*x)^2)) +
          1/(4*k)*(2*k+4*k^2*x/sqrt(1+(2*k*x)^2))/(sqrt(1+(2*k*x)^2)+2*k*x)
  # pure elliptical exception (no camber)
  if k == 0
    s(x) = x
    ds(x) = 0
  end

  s_c2 = s(c/2)
  #print(s_c2,"\n")
  s_max = 2*s_c2

  # calculate camber projection point and the the new shape cooinate based on elliptical
  for i = 1:num_p+1

    x = xy[i,1]
    s_real(a) = x/c*s_max-s_c2 - s(a) # function for newton iteration
    ds_real(a) = -1*ds(a) #
    x_c = c/2 + Newton_method(s_real,ds_real,x-c/2) # camber x position
    y_c = y_para(x_c)
    dy_c = dy_para(x_c)

    t_2 = xy[i,2]   # thickness/2
    L = sqrt(1 + dy_c^2)
    xy[i,1] = x_c - dy_c/L*t_2
    xy[i,2] = y_c + 1/L*t_2
    xy_c[i,1] = x_c
    xy_c[i,2] = y_c

  end

  # Outpit file and distribute to upper surface and lower surface

  xyz_u[1:num_up,1:2] = xy[1:num_up,:]
  xyz_l[1:num_lp-1,1:2] = xy[num_up:num_p,:]
  xyz_l[num_lp,1:2] = xy[1,:]
  # write the file output
  # out_type = 1:xfoil data; = 2: 3 section data
  if out_type == 1

     fid = open("ellip_t$(max_t*100)c$(camber*100).xy","w")
     println(fid,"ellip_t$(max_t*100)c$(camber*100)")
     writedlm(fid,xy)
     writedlm(fid,xy[1,:]')
     close(fid)

  elseif out_type == 2

     fid = open("ellip_t$(max_t*100)c$(camber*100)_us.dat","w")
     writedlm(fid,num_up)
     writedlm(fid,xyz_u)
     close(fid)

     fid = open("ellip_t$(max_t*100)c$(camber*100)_ls.dat","w")
     writedlm(fid,num_lp)
     writedlm(fid,xyz_l)
     close(fid)

  end

  # plot the shape and camber profile
  pl = plot(xy[:,1],xy[:,2],aspect_ratio = 1,
      title = "ellip_t$(max_t*100)c$(camber*100)",labels ="shape")
  pl = plot!(xy_c[:,1],xy_c[:,2],labels = "camber")
  display(pl)

end

function refinement(r::Float64,P1::Int64,P2::Int64,theta::Vector{Float64})

     n = abs(P2-P1)  #P2 P1 is the data point index in series, n will be section number
                     # P2 is the coarse side P1 is the fine side
     b = zeros(Float64,n) #b is tje new section length after refinement

     T = abs(theta[P2]-theta[P1]) # total length of the refinement section
     dir =convert(Int64,(P2-P1)/abs(P2-P1)) # dir is used to define whether distance should
                             # be added forward the point or subtracted backward

     for i = 1:n
       b[i] = T*(r-1)/(r^n-1)*r^(i-1)
       theta[P1+dir*i] = theta[P1+dir*(i-1)] + dir*b[i]
     end
  return theta

end

function refinement_iteration(x::Function,y::Function,num_p::Int64,r_p::Float64)

  # r_p    = 0.25 # refinement portion : should be 0 ~ 25%
  #r_ratio = 1.01 # refinement gepmetric seroies ratio, can not be 1.0
                  # since it will break the geometric seried

  # initilize equal space input for doorcinate
  theta = zeros(Float64,num_p+1)

  for i  = 1:num_p+1
     theta[i] = 2*pi/num_p*(i-1)
  end

  # iteration to find r_ratio for the segment near the edge is small enough
  xy_edge = zeros(Float64,2,2) # initialize xy data near edge

  for iter = 1 : 100
    r_ratio = 1.01 + (iter-1)*0.01
    theta = refinement(r_ratio,1,round(Int64,num_p*r_p),theta)

    for i  = 1:2
    xy_edge[i,1] = x(theta[i])
    xy_edge[i,2] = y(theta[i])
    end

    e_l = norm(xy_edge[2,:]-xy_edge[1,:]) #near edge section length
    if e_l <= 1e-5
      print("near edge length = ",e_l,"\n")
      print("refinement ratio=",r_ratio,"\n")
      print("refinement succeed \n")
      return r_ratio
    end
  end
  @warn("refinement failed, can't get sufficient
         refinement_ratio under the refinement portion and total point setting")

end

function Newton_method(f::Function,f0::Function,x0::Float64,
  tol::Float64 = 1e-5,maxiter::Integer=200, eps::Float64=1e-10,m::Float64 = 0.7)

     for i = 1:maxiter
       f0_value = f0(x0)
       if abs(f0_value) < eps
          # @warn("first derivative is zero! \n")
          return x0
       end
       f_value = f(x0)
       x_n1 = x0 - m*f_value/f0_value
       if abs((x_n1-x0)/x0) < tol
          print("sol =",x_n1,"\n")
          return x_n1
       end
       x0 = x_n1
     end
    print("sol =",x_n1,"\n")
    @error("max iteration reached \n")
end

# NACA shape generator
#=
function NACA_4()


  num_tp = round(Int64,num_p*0.03) # TE point number (3% ofoverall point)
  num_up = round(Int64,(num_p-num_tp+2)/2) # upper surface point number
  num_lp = round(Int64,num_p-num_tp+1-num_up+2) # lower surface point number

  xy  = zeros(Float64,num_p,2)    # xy data set (for xfoil))
  xy_c = zeros(Float64,num_p,2)   # xy data for camber line
  xyz_l = zeros(Float64,num_lp,3) # xyz data set for lower surface
  xyz_u = zeros(Float64,num_up,3) # xyz data set for upper surface
  xyz_t = zeros(Float64,num_tp,3) # xyz data set for trailing edge section


  # refinement iteration
  x_naca(a) = c/2*(1+cos(a)) # elliptical x function (a= theta here)
  y_naca(a) = max_t/2*sin(a) # elliptical y function (a= theta here)
  r_p = 0.25  # refinement portion, should be between 0~25%
  r_ratio = refinement_iteration(x_ellip,y_ellip,num_p,r_p)

   # initialize the theta
   for i  = 1:num_p+1
      theta[i] = 2*pi/num_p*(i-1)
   end
   #  do the actual refinement
   theta = refinement(r_ratio,1,round(Int64,num_p*r_p),theta)
   theta = refinement(r_ratio,round(Int64,num_p*0.5),round(Int64,num_p*(0.5-r_p)),theta)
   theta = refinement(r_ratio,round(Int64,num_p*0.5),round(Int64,num_p*(0.5+r_p)),theta)
   theta = refinement(r_ratio,num_p+1,round(Int64,num_p*(1-r_p)),theta)


   # distribute the data to xy_l,xy_u and xy_t section
   p_start = num_tp-round(Int8,num_tp/2)
   p_end = p_start + num_up-1
   xyz_u[1:num_up,1:2] = xy[p_start:p_end,:]

   p_start = p_end
   p_end = p_start + num_lp-1
   xyz_l[1:num_up,1:2] = xy[p_start:p_end,:]

   p_start = p_end
   p_end = p_start + num_tp-1
   for i = p_start:p_end
       j = i + 1 - p_start
       z = mod(i,num_p)
       if z == 0
          z = num_p
       end
       xyz_t[j,1:2] = xy[z,:]
   end

   # write the file output
   # out_type = 1:xfoil data; = 2: 3 section data
   if out_type == 1

      fid = open("ellip_t$(max_t*100)c$(camber*100).xy","w")
      println(fid,"ellip_t$(max_t*100)c$(camber*100)")
      writedlm(fid,xy)
      writedlm(fid,xy[1,:]')
      close(fid)

   elseif out_type == 2

      fid = open("ellip_t$(max_t*100)c$(camber*100)_te*.dat","w")
      writedlm(fid,num_tp)
      writedlm(fid,xyz_t)
      close(fid)

      fid = open("ellip_t$(max_t*100)c$(camber*100)_us*.dat","w")
      writedlm(fid,num_up)
      writedlm(fid,xyz_u)
      close(fid)

      fid = open("ellip_t$(max_t*100)c$(camber*100)_ls*.dat","w")
      writedlm(fid,num_lp)
      writedlm(fid,xyz_l)
      close(fid)

   end
end
=#
main()
