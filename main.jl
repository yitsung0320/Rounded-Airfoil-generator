# Rounded Edge Airfoil generator
# Developed by Yi Tsung Lee, NCSU AE
# The generated file will have the LE coordinate near (0,0)
# The generated file will have the TE coordinate at (c,0)
# change the

using Plots
using DelimitedFiles
gr()

function main()

   # =============== control variable section ================
   # out_type = 1:xfoil data form; = 2:3_section data
   out_type = 2

   # variable definition
   c::Float64      = 1.0 # chord length
   max_t::Float64  = 0.12 # max thickness (between 0~1.0)
   max_tp::Float64 = 0.0 # max thickness position(in x/c,should be 0~1.0)
   num_p::Int64    = 200 # number of points
   camber::Float64 = 0.06 # camber (in y/c )
   # =============== control variable section ================

   #initilize
   num_tp = round(Int8,num_p*0.03) # TE point number (3% ofoverall point)
   num_up = round(Int8,(num_p-num_tp+2)/2) # upper surface point number
   num_lp = round(Int8,num_p-num_tp+1-num_up+2) # lower surface point number

   xy  = Array{Float64,2}(undef,num_p,2)   # xy data set
   xy_l = Array{Float64,2}(undef,num_lp,2) # xy data set for lower surface
   xy_u = Array{Float64,2}(undef,num_up,2) # xy data set for upper surface
   xy_t = Array{Float64,2}(undef,num_tp,2) # xy data set for trailing edge section

   # start of the coordinate generation

   xy = ellipse(num_p,c,max_t,camber,xy)
   pl = plot(xy[:,1],xy[:,2],aspect_ratio = 1,
        title = "ellip_t$(max_t*100)c$(camber*100)",labels = "shape")
   display(pl)

   # distribute the data to xy_l,xy_u and xy_t section
   p_start = num_tp-round(Int8,num_tp/2)
   p_end = p_start + num_up-1
   xy_u = xy[p_start:p_end,:]

   p_start = p_end
   p_end = p_start + num_lp-1
   xy_l = xy[p_start:p_end,:]

   p_start = p_end
   p_end = p_start + num_tp-1
   for i = p_start:p_end
       j = i + 1 - p_start
       z = mod(i,num_p)
       if z == 0
          z = num_p
       end
       xy_t[j,:] = xy[z,:]
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

      fid = open("ellip_t$(max_t*100)c$(camber*100)_te.xy","w")
      println(fid,"ellip_t$(max_t*100)c$(camber*100)_te")
      writedlm(fid,xy_t)
      close(fid)

      fid = open("ellip_t$(max_t*100)c$(camber*100)_us.xy","w")
      println(fid,"ellip_t$(max_t*100)c$(camber*100)_us")
      writedlm(fid,xy_u)
      close(fid)

      fid = open("ellip_t$(max_t*100)c$(camber*100)_ls.xy","w")
      println(fid,"ellip_t$(max_t*100)c$(camber*100)_ls")
      writedlm(fid,xy_l)
      close(fid)

   end

end


function ellipse(num_p::Int64,c::Float64,max_t::Float64,
                     camber::Float64,xy::Array{Float64,2})

 # ellipse based coordinate
  for i  = 1:num_p
     theta = 2*pi/num_p*(i-1)
     xy[i,1] = c/2*(1+cos(theta))
     xy[i,2] = max_t/2*sin(theta)
  end
 # apply parabolic camberline
  for i = 1:num_p
    dy = -4*camber/c*xy[i,1]*(xy[i,1]-c)
    xy[i,2] = xy[i,2] + dy
  end

 return xy
end

# NACA shape generator
function NACA_4()
end

main()
