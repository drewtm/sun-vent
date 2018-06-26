#math is standard with python
from math import pi, tan, sin, cos, atan, sqrt
#ezdxf is from https://github.com/mozman/ezdxf
from ezdxf.r12writer import r12writer
#svgwrite is from https://github.com/mozman/svgwrite
import svgwrite
from svgwrite import mm

#A series of concentric, conical mirrors arranged in a disc can focus light like a lens,
#with the focal area on the opposite side of the disc from the light source.
#This script calculates the geometry of a series of cones with specified geometric parameters
#and then generates a flat pattern of annular arcs which can be bent into those cones.
#The calculations assume a far away light source with parallel incoming rays, like the sun.

#An assumption is made that these are front-face mirrors, and the finite-thickness edge of
#each cone is normal to the axis of the array. For thick cone material, rear-face mirrors,
#wide arrays, and/or short focal lengths, this approximation may be poor.

#the efficiency and focusing power calculations assume that light from each cone is spread
#evenly over its whole focal disc, which untrue. However, the real distribution
#of incident light flux in the focal disc is difficult to represent and probably would
#not yield accurate results anyway, given the probable construction quality of these arrays

#maximum height of any of the conical tubes
hmax = 50.0

#limiting diameter of the widest cone (the finished array may be a bit smaller than this)
D_max = 200.0
r_max = D_max/2 #limiting radius for the array

#distance from the mid-plane of the array to the focal plane
f = 200.0

#maximum size of the focused light area
D_f_max = 10.0
r_f_max = D_f_max/2 #max radius of the focused light area

#thickness of the sheet from which the mirrors are cut
thk = 0.4

#kerf allowance for the material that is removed around the part when flat pattern is cut from the sheet material
kerf = 1.0

#how accurately can you align the array toward the light source?
#this will define the angle of the steepest middle cone, which is
#the steepest face inclination measured from the midplane of the array
angle_tolerance = 2.0 #degrees

#----------ok, lets do the math

#innermost face inclination in radians
q1 = (90-angle_tolerance)*pi/180

#the radius of the innermost cone is set to direct the
#light that hit the cone mid-face directly onto the focal point
r1 = f/(tan(2*q1-pi/2))

#top and bottom radii of the innermost cone are based on
#extrapolating up and down to the constrained height of the cone
r1_bot = r1 - hmax/2/tan(q1)
r1_top = r1 + hmax/2/tan(q1)

#face height of the innermost cone, slightly longer than the height of the array
#because the face is tilted
L1 = hmax/sin(q1)

#radius of the light blob projected by the innermost cone
r_foc = L1/2*(sin(q1)*r1/f-cos(q1))

#initialize lists for the cone geometries
r_bot = [r1_bot] #bottom/small radius
r_top = [r1_top] #top/large radius
q = [q1] #face angle
L = [L1] #face length
H = [hmax] #axial height
#radius on the focal plane in which the light from this cone is projected
r_out = [r_foc]
#incoming light area that each cone receives, used for efficiency calculations
area_in = [(r1_top*r1_top - r1_bot*r1_bot)*pi]
#area on the focal plane in which the light from this cone is projected
area_out = [r_foc*r_foc*pi]

#now work outward through the cones, stopping at the maximum radius.
#each cone is sized so that the light that comes from its top edge will
#barely skim under the bottom of the previous cone to the focal point,
#and that the light bouncing off the bottom edge won't fall outside the
#maximum focal area radius
while True:
	#angle between the midplane of the array and the light skimming
	#past the bottom of the previous ring to the focal point
	wn = atan((f - hmax/2) / (r_bot[-1] + thk))
	
	#face inclination for this cone
	qn = pi/4 + wn/2
	
	#the in-plane projected radius of the cone that completely fills the focal disc
	#with reflected light, and light from the top of the cone just skirts the previous cone
	delta_rn = (r_f_max)*tan(wn)/(tan(qn)-tan(wn))
	
	#the height of the cone with delta_rn
	cn = delta_rn*tan(qn)
	
	#if the cone with height cn would be too tall, then we need to shorten it
	#to a height of h and shrink its radius until its light just skirts the previous cone
	if cn > hmax:
		#height of this cone is now defined as hmax
		Hn = hmax
	
		#top radius of the cone based on the max height
		rn_top = r_bot[-1] + thk + hmax/tan(wn)
		
		#bottom radius of the cone:
		rn_bot = rn_top - hmax/tan(qn)
		
		#focal radius of this cone
		rf = rn_bot - (f - hmax/2)/tan(wn)
		
	#otherwise, we'll leave this cone shorter than hmax
	else:
		#height of this cone
		Hn = cn
		
		#bottom radius that will reflect light to the edge of the focus disc
		rn_bot = r_bot[-1] + thk + r_f_max
		
		#top radius of the cone based on the max focal disc radius
		rn_top = rn_bot + delta_rn
		
		#focal radius of this cone is still the maximum allowed
		rf = r_f_max
	
	#check to see if this cone exceeds the maximum allowed diameter
	if rn_top > r_max :
		break
	
	#face width of this cone
	Ln = (rn_top-rn_bot)/cos(qn)
	
	#axially projected area of incoming light that this cone consumes
	an_in = (rn_top*rn_top - rn_bot*rn_bot)*pi

	#area on the focal plane that this cone projects into
	an_out = rf*rf*pi
	
	#now add this ring to the output lists
	r_top.append(rn_top)
	r_bot.append(rn_bot)
	q.append(qn)
	L.append(Ln)
	H.append(Hn)
	r_out.append(rf)
	area_in.append(an_in)
	area_out.append(an_out)

#reverse the order of the lists, starting from the outermost cone
r_top.reverse()
r_bot.reverse()
q.reverse()
L.reverse()
H.reverse()
r_out.reverse()
area_in.reverse()
area_out.reverse()

#---------now lets calculate some statistics about this design

#the total axially projected area of light captured by cones
a_captured = sum(area_in)
#don't forget that the light is also shining directly on the focal area, through the middle!
#but lets only count the light falling on already brightened areas for efficiency calculation
a_captured += area_out[0]
#the total footprint of the cone array
a_total = r_top[0]*r_top[0]*pi
#what percentage of light hit a cone face (geometric efficiency)
area_efficiency = a_captured/a_total

#initialze some arrays
cumulative_efficiency = [] #the efficiency of the array if this were the innermost ring
cumulative_area = [] #the projected area of the array rings if this were the innermost ring
#find the overall efficiency added by each ring, starting from the outside
#to help identify how useful the small inner rings are
for i, a in enumerate(area_in):
	cumulative_area.append(a+sum(area_in[0:i]))
	cumulative_efficiency.append(cumulative_area[i]/a_total)

#find the average flux factor over the entire focused area
#this would be multiplied with the flux of the incident light to find
#the average flux in the focal area, assuming perfectly reflective surfaces
avg_flux_factor = 0
for a_in, a_out in zip(area_in, area_out):
	primary = a_in/area_out[0] #how much the light is concentrated by this ring
	avg_flux_factor += a_in*primary #sum is weighted by the area of each ring
avg_flux_factor /= (a_captured-area_out[0]) #divide the total projected area back out of the sum
avg_flux_factor += 1 # add in the unconcentrated light shining through the opening of the middle cone

#calculate the flat pattern dimensions
r_i_flat = [] #inner radius of flattened arc
r_o_flat = [] #outer radius of flattened arc
ang_flat = [] #included angle of flattened arc
for rt, rb, dr in zip(r_top, r_bot, L):
	#when the arcs are bent into cones, the inner and outer arc length become the top and bottom
	#circumference, but this is measured in the neutral axis halfway through the material thickness
	#while the calculated optical surface is the inner surface of the cone. so we have to add half the
	#material thickness to the calculated radii to make the flat pattern the right size.
	circ_i = (rb+thk/2)*2*pi #inner circumference becomes the inner arc length
	circ_o = (rt+thk/2)*2*pi #outer circumference becomes the outer arc length
	ang_flat.append((circ_o - circ_i)/dr) #the face width does not change from flat arc to bent cone
	r_i_flat.append(circ_i/ang_flat[-1]) #flat arc radius can be calculated from the included angle and the arc length
	r_o_flat.append(circ_o/ang_flat[-1]) #same for outer

#-----------now lets make some useful outputs!
	
#code to write a dxf flat pattern cut file
with r12writer("arcs.dxf") as dxf:
	r_prev = 0 #initialize a helper variable to remember the radius of the next largest ring
	left_end = 0 #initialize a helper variable to remember the leftmost point of the series of arcs in the layout
	for ri, ro, a in zip(r_i_flat, r_o_flat, ang_flat): #again, this is stepping through the rings from largest to smallest
		#first we do a little math to space the rings nicely together in the flat pattern
		hx = ro*sin(a/2) #the height of the corner where this ring comes closest to the next largest ring on the flat pattern
		if r_prev == 0: #don't try to take the square root of a negative number!
			xx = 0
		else:
			#a horizontal distance from the previous arc's center to the outer corner overlap point
			x_prev = sqrt(r_prev*r_prev - hx*hx)
			#the new horizontal position of this arc's center
			xx += (x_prev-ro*cos(a/2)-kerf)
		dxf.add_arc((xx, 0), radius=ri, start=-a*90/pi, end=a*90/pi) #draw the inner arc
		dxf.add_line((xx+ri*cos(a/2), ri*sin(a/2)), (xx+ro*cos(a/2), ro*sin(a/2))) #draw the top butt-joint edge
		dxf.add_arc((xx, 0), radius=ro, start=-a*90/pi, end=a*90/pi) #draw the outer arc
		dxf.add_line((xx+ro*cos(a/2), -ro*sin(a/2)), (xx+ri*cos(a/2), -ri*sin(a/2))) #draw the bottom butt-joint edge
		left_end = xx+ri*cos(a/2)
		r_prev = ri #remember how big this inner arc was for the next time through the loop
		
	#now draw a pair of straight, slotted strips that can be used to hold all the rings concentric
	ds = max(thk*5.0, hmax/8.0) #depth of slots is 1/8 of the array height or 5x the thickness of the material
	ws = float(max(thk, kerf)) #don't make the slots thinner than the allowable kerf, such that they're uncuttable
	dh = ds + r_top[0]/10.0 #height of strip below the slots is 1/20 of the length of the strip
	y_slotted_edge = r_o_flat[0]*sin(ang_flat[0]/2) - dh
	x_offset = left_end+r_top[0]
	pt1 = (x_offset+r_top[0], y_slotted_edge)
	pt2 = (pt1[0], pt1[1]+dh)
	dxf.add_line(pt1, pt2) #right vertical edge of strip
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
	pt1 = pt2
	pt2 = (x_offset+ws/2, pt1[1])
	dxf.add_line(pt1, pt2) #over to center slot
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
	pt1 = pt2
	pt2 = (pt1[0], pt1[1]-dh/2)
	dxf.add_line(pt1, pt2) #down into slot
	pt1 = pt2
	pt2 = (pt1[0]-ws, pt1[1])
	dxf.add_line(pt1, pt2) #bottom of slot
	dxf.add_line((pt1[0], pt1[1]-dh/2-kerf), (pt2[0], pt2[1]-dh/2-kerf)) #second copy but no slot
	pt1 = pt2
	pt2 = (pt1[0], pt1[1]+dh/2)
	dxf.add_line(pt1, pt2)
	pt1 = pt2
	pt2 = (x_offset-r_top[0], pt1[1])
	dxf.add_line(pt1, pt2)
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
	pt1 = pt2
	pt2 = (pt1[0], y_slotted_edge)
	dxf.add_line(pt1, pt2)
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
	pt1 = pt2
	for rb, a in zip(r_bot, q):
		pt2 = (x_offset-rb-ws-ds/tan(a), pt1[1])
		dxf.add_line(pt1, pt2) #line from previous slot to this slot along slotted edge
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]), (-pt2[0]+2*x_offset, pt2[1])) #and its mirror
		dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]-dh-kerf), (-pt2[0]+2*x_offset, pt2[1]-dh-kerf)) #mirror for 2nd copy
		pt1 = pt2
		#depth of each slot is increased to allow for the bit of cone material thickness that
		#protrudes axially past the optical edge
		pt2 = (pt1[0]+ds/tan(a), pt1[1]+ds+thk*cos(a))
		dxf.add_line(pt1, pt2) #line from slotted edge down into slot
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]), (-pt2[0]+2*x_offset, pt2[1])) #and its mirror
		dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]-dh-kerf), (-pt2[0]+2*x_offset, pt2[1]-dh-kerf)) #mirror for 2nd copy
		pt1 = pt2
		pt2 = (pt1[0]+ws, pt1[1])
		dxf.add_line(pt1, pt2) #line in base of slot
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]), (-pt2[0]+2*x_offset, pt2[1])) #and its mirror
		dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]-dh-kerf), (-pt2[0]+2*x_offset, pt2[1]-dh-kerf)) #mirror for 2nd copy
		pt1 = pt2
		pt2 = (pt1[0], y_slotted_edge)
		dxf.add_line(pt1, pt2) #line coming straight back up out of the slot
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]), (-pt2[0]+2*x_offset, pt2[1])) #and its mirror
		dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
		dxf.add_line((-pt1[0]+2*x_offset, pt1[1]-dh-kerf), (-pt2[0]+2*x_offset, pt2[1]-dh-kerf)) #mirror for 2nd copy
		pt1 = pt2
	pt2 = (x_offset-ws/2, y_slotted_edge)
	dxf.add_line(pt1, pt2)
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
	pt1 = pt2
	pt2 = (pt1[0], pt1[1]+dh/2)
	dxf.add_line(pt1, (pt2[0], pt2[1]-dh/2)) #no slot on this one
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy but with slot
	pt1 = pt2
	pt2 = (pt1[0]+ws, pt1[1])
	dxf.add_line((pt1[0], pt1[1]-dh/2), (pt2[0], pt2[1]-dh/2)) #no slot on this one
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy but with slot
	pt1 = pt2
	pt2 = (pt1[0], y_slotted_edge)
	dxf.add_line((pt1[0], pt1[1]-dh/2), pt2)
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy
	pt1 = pt2
	pt2 = (x_offset+rb, y_slotted_edge)
	dxf.add_line(pt1, pt2)
	dxf.add_line((pt1[0], pt1[1]-dh-kerf), (pt2[0], pt2[1]-dh-kerf)) #second copy

#this creates a preview of what the finished array should look like and its properties
im = svgwrite.Drawing("cones.svg", size=((D_max+20)*mm, (r_max+hmax+20)*mm), debug=True)
im.viewbox(-(D_max+20)*2.5, -(hmax+10)*5, (D_max+20)*5, (r_max+hmax+20)*5)
#draw the rings, with the black stroke representing the material thickness
for ro, ri, h in zip(r_top, r_bot, H):
	im.add(im.circle(center=(0,0),r=(ro+thk/2)*mm,fill='gray', stroke='black', stroke_width=thk*mm))
	im.add(im.circle(center=(0,0),r=(ri*mm), fill='white'))
#hide the upper half to draw the cross-section
im.add(im.rect(insert=(-r_max*mm, -r_max*mm), size=(D_max*mm, r_max*mm), fill='white'))
#draw a cross-section of the rings in the top half of the image
for ro, ri, h, rf in zip(r_top, r_bot, H, r_out):
	im.add(im.line(start=((ri+thk/2)*mm, 0), end=((ro+thk/2)*mm, -h*mm), stroke='black', stroke_width=thk*mm))
	im.add(im.line(start=(-(ri+thk/2)*mm, 0), end=(-(ro+thk/2)*mm, -h*mm), stroke='black', stroke_width=thk*mm))
	#include a shaded preview of the quality of the focal area
	im.add(im.circle(center=(0, 0), r=rf*mm, fill='red', opacity=(1.0/len(H))))
#one more layer to represent the unfocused light
im.add(im.circle(center=(0, 0), r=r_bot[-1]*mm, fill='red', opacity=(1.0/len(H))))
#hide the top half of the shaded focus preview
im.add(im.rect(insert=(-r_f_max*mm,-r_f_max*mm), size=(D_f_max*mm, r_f_max*mm), fill='white'))
#draw a height graph of sorts showing how the light is distributed in the focal area
for rf, hbar, tbar in zip(r_out, cumulative_area, area_in):
	#the thickness of each of these bars corresponds to the area of the ring that added the light
	im.add(im.line(start=(-rf*mm, (-hbar+area_in[0]/2)*hmax/a_captured*mm), end=(rf*mm, (-hbar+area_in[0]/2)*hmax/a_captured*mm), stroke='red', stroke_width=tbar*hmax/a_captured*mm))
#now add some useful statistics and performance numbers about the array
im.add(im.text('{:d} rings'.format(len(H)), (0, -(hmax+13.5)*mm), font_size=3*mm, text_anchor='middle'))
im.add(im.text('average focused brightness is {:.0f}x'.format(avg_flux_factor), (0, -(hmax+6.5)*mm), font_size=3*mm, text_anchor='middle'))
im.add(im.text('{:.0f}% of light is captured'.format(area_efficiency*100), (0, -(hmax+10)*mm), font_size=3*mm, text_anchor='middle'))
im.save()