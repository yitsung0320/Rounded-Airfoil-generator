# Rounded Edge Airfoil generator
# Developed by Yi Tsung Lee, NCSU AE
# The generated file will have the LE coordinate near (0,0)
# The generated file will have the TE coordinate at (c,0)
# outputfile will be in xfoil format(out_type =1)
# outputfile include num point in the section and 3-d coordinate (out_type =2)

# Remenber to check the working directory first for file saving

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

   xy  = zeros(Float64,num_p,2)   # xy data set (for xfoil))
   xyz_l = zeros(Float64,num_lp,3) # xyz data set for lower surface
   xyz_u = zeros(Float64,num_up,3) # xyz data set for upper surface
   xyz_t = zeros(Float64,num_tp,3) # xyz data set for trailing edge section

   # start of the coordinate generation

   xy = ellipse(num_p,c,max_t,camber,xy)
   pl = plot(xy[:,1],xy[:,2],aspect_ratio = 1,
        title = "ellip_t$(max_t*100)c$(camber*100)",labels = "shape")
   display(pl)

   # distribute the data to xy_l,xy_u and xy_t section
   p_start = num_tp-round(Int8,num_tp/2)
   p_end = p_start + num_up-1
   xyz_u[1:num_up,1:2] = xy[p_start:p_end,:]

   p_start = p_end
   p_end = p_start + num_lp-1
   xyz_l = xy[p_start:p_end,:]

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

      fid = open("ellip_t$(max_t*100)c$(camber*100)_te.xy","w")
      writedlm(fid,num_tp)
      writedlm(fid,xyz_t)
      close(fid)

      fid = open("ellip_t$(max_t*100)c$(camber*100)_us.xy","w")
      writedlm(fid,num_up)
      writedlm(fid,xyz_u)
      close(fid)

      fid = open("ellip_t$(max_t*100)c$(camber*100)_ls.xy","w")
      writedlm(fid,num_lp)
      writedlm(fid,xyz_l)
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
