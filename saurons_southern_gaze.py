import numpy as np
import sys
fo = "saurons_southern_gaze.svg"
fi = "southern_begin.txt"
with open(fi, "r") as tf:
  str_begin = tf.read()
fi = "southern_end.txt"
with open(fi, "r") as tf:
  str_end = tf.read()

cx = 1361.26
cy = 747.37-(747.37-715.59)
cxl = 311.26001
cxr = 2411.26
cxof = cxr - cxl
# add pupil:
cxo = 2411.26 - 1361.26
cxoo = cxo #+ ((1120 - cxo) % pra[1])
pr = np.linspace(0,650, int(np.floor(156+1)))
sr = pr[2]
pra = np.linspace(0,650, int(np.floor(156/2)+1))
pre = np.arange(cxoo, 1123, pra[1]/2)
R = pre[-1]
L = cxof/2
H = np.sqrt(R**2 - (cxof/2)**2)
lwidth = 0.1

annotate_index = sys.argv[1]
if (sys.argv[1] == "-i"):
  annotate_index = 0
else:
  annotate_index = sys.argv[1]
  
print("rows in: {}".format(annotate_index))
print("rows in: {}".format(annotate_index), file=open("southern_progress.txt", "w"))
with open(fo, "w") as text_file:

  print(str_begin, file=text_file)

  str_next = "<path\n     style=\"fill:none;stroke:#ffffff;stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n     d=\"m {},{} {},{} z\"\n     id=\"path{}\"\n     inkscape:connector-curvature=\"0\" />".format(lwidth*4, cx, cy-H, 0, 2*H, 124412)
  print(str_next, file=text_file)

  i = 5183
  # add stitch lines to eye segments
  
  jjj = 100001
  print("pre : {}".format(pre))
  for r in pre[1:]:
    j = int(np.floor(r/sr)+1)
    jj = int(np.floor(r/(sr/2) + 1))
    if (jjj - 100000 < int(annotate_index)):
      print("Printing Red Ring {}".format(j))
      colore="#ff0000"
    else:
      colore="#ffffff"

    print("j: {}, r/sr: {}".format(j, r/sr))
    if True: #np.abs(j - r/sr) > 0.7:
      print("j: {} should get a line!".format(j, r/sr))
      dangle = (2*np.pi)/(j*10+1)
      angles = np.arange(0, np.arccos(L/r), dangle)
      #angles = np.linspace(0, 
      k = 0
      for angle in angles:
        for sign in [-1, -1]:
          angle = sign*angle

          mod_factor = 1
          #if (k % 6 == 0):
          #  mod_factor = 2
          if (k % j == 0):
            mod_factor = 4
          if (k == 0):
            mod_factor = 8
          if (k % 10 == 0):
            mod_factor = 4
          # check that the distance from each eye bounding circle center to end point is less than R
          #cy = 747.37-(747.37-715.59)
          cxl = 311.26001
          cxr = 2411.26

          xol = cxl+(jj-1)*sr*np.cos(angle)/2
          yol = cy+(jj-1)*sr*np.sin(angle)/2
          xpl = -sr*np.cos(angle)/2
          ypl = -sr*np.sin(angle)/2

          xor = cxr-(jj-1)*sr*np.cos(angle)/2
          yor = cy+(jj-1)*sr*np.sin(angle)/2
          xpr = sr*np.cos(angle)/2
          ypr = -sr*np.sin(angle)/2

          # oouter-end
          dl = np.sqrt((xol+xpl-cxl)**2 + (yol+ypl-cy)**2)
          dr = np.sqrt((xor+xpr-cxr)**2 + (yor+ypr-cy)**2)
          # inner-end
          dli = np.sqrt((xol-cxl)**2 + (yol-cy)**2)
          dri = np.sqrt((xor-cxr)**2 + (yor-cy)**2)
          # check if end points are left (resp right) of vertical
          if xor < cx:
            # if left end point is, but left starting point is not, then chop
            if (xor+xpr > cx):
              # compute intersection with vertical
              dr = r - L/np.cos(angle)
              xpr = dr*np.cos(angle)
              ypr = -dr*np.sin(angle)
              # chop
          if xol > cx:
            if (xol+xpl < cx):
              # compute intersection with vertical
              dr = r - L/np.cos(angle)
              xpl = -dr*np.cos(angle)
              ypl = -dr*np.sin(angle)
              # chop

            str_next = "<path\n     style=\"fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n     d=\"m {},{} {},{} z\"\n     id=\"pathl{}\"\n     inkscape:connector-curvature=\"0\" />".format(colore, lwidth*mod_factor, xol, yol, xpl, ypl, i)
            print(str_next, file=text_file)

            str_next = "<path\n     style=\"fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n     d=\"m {},{} {},{} z\"\n     id=\"pathr{}\"\n     inkscape:connector-curvature=\"0\" />".format(colore, lwidth*mod_factor, xor, yor, xpr, ypr, i)
            print(str_next, file=text_file)
        

        i+=1
        k+=1

        
    print("generating eye segments of radius r: {}".format(r))
    mod_factor = 1
    if (jjj % 2 == 0):
      mod_factor = 2.0
    if (jjj % 10 == 0):
      mod_factor = 4.0
    if (jjj % 20 == 0):
      mod_factor = 8.0

    # left side
    offset = np.arccos(cxof/2/r)
    la_b = np.pi+offset
    la_e = np.pi-offset
    ra_b = offset
    ra_e = 2*np.pi-offset
    # right side
    print("generating eye segments of radius r: {}".format(r))
    str_before_l =" <path\n      sodipodi:cx=\"311.26001\"\n sodipodi:cy=\"715.59003\"\n sodipodi:rx=\"{}\"\n sodipodi:ry=\"{}\"\n sodipodi:end=\"{}\"\n sodipodi:start=\"{}\"\n sodipodi:open=\"true\"\n sodipodi:type=\"arc\"\n d=\"m 1350.0359,296.84816 a {},{} 0 0 1 -7.1876,854.89044\"\n id=\"path31910\"\n style=\"fill:none;fill-opacity:1;stroke:{};stroke-width:{};stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1\" />".format(r, r, ra_b, ra_e, r, r, colore, mod_factor*0.1)

    str_before_r = "<path \n style=\"fill:none;fill-opacity:1;stroke:{};stroke-width:{};stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1\"\n id=\"path31893\"\n d=\"M 1376.0413,1143.0503 A {},{} 0 0 1 1362.4285,322.71281\"\n sodipodi:type=\"arc\"\n sodipodi:open=\"true\"\n sodipodi:start=\"{}\"\n sodipodi:end=\"{}\"\n sodipodi:ry=\"{}\"\n sodipodi:rx=\"{}\"\n sodipodi:cy=\"715.59003\"\n sodipodi:cx=\"2411.26\" /> ".format(colore, mod_factor*0.1, r, r, la_e, la_b, r, r)

    print(str_before_l, file=text_file)
    print(str_before_r, file=text_file)    
    i+=1
    jjj+=1


  j = 1
  for r in pr[1:]:
    #print(r)
    mod_factor = 1
    if (j % 2 == 0):
      mod_factor = 2.0
    if (j % 10 == 0):
      mod_factor = 4.0
    if (j % 20 == 0):
      mod_factor = 8.0

    if (j < int(annotate_index)):
      print("Printing Red Ring {}".format(j))
      colore="#ff0000"
    else:
      colore="#ffffff"


    if (r < H) & (r > R - cxof/2):
      offset = np.pi - np.arccos((r**2+(cxof/2)**2 - R**2)/(r*cxof))
      la_b = np.pi+offset
      la_e = np.pi-offset
      ra_b = offset
      ra_e = 2*np.pi-offset
      # right side
      str_before_l =" <path\n      sodipodi:cx=\"1361.26\"\n sodipodi:cy=\"715.59003\"\n sodipodi:rx=\"{}\"\n sodipodi:ry=\"{}\"\n sodipodi:end=\"{}\"\n sodipodi:start=\"{}\"\n sodipodi:open=\"true\"\n sodipodi:type=\"arc\"\n d=\"m 1350.0359,296.84816 a {},{} 0 0 1 -7.1876,854.89044\"\n id=\"path{}\"\n style=\"fill:none;fill-opacity:1;stroke:{};stroke-width:{};stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1\" />".format(r, r, ra_b, ra_e, r, r, j, colore, mod_factor*0.1)

      str_before_r = "<path \n style=\"fill:none;fill-opacity:1;stroke:{};stroke-width:{};stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1\"\n id=\"path{}\"\n d=\"M 1376.0413,1143.0503 A {},{} 0 0 1 1362.4285,322.71281\"\n sodipodi:type=\"arc\"\n sodipodi:open=\"true\"\n sodipodi:start=\"{}\"\n sodipodi:end=\"{}\"\n sodipodi:ry=\"{}\"\n sodipodi:rx=\"{}\"\n sodipodi:cy=\"715.59003\"\n sodipodi:cx=\"1361.26\" /> ".format(colore, mod_factor*0.1, j, r, r, la_e, la_b, r, r)
      print(str_before_l, file=text_file)
      print(str_before_r, file=text_file)    

    elif (r > H):
      str_before = "  <circle \n     style=\"fill:none;stroke:{};stroke-width:{};stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1\" \n     id=\"circle{}\"\n     cx=\"{}\"\n     cy=\"{}\"\n     r=\"{}\" />".format(colore, mod_factor*0.1, i, cx, cy, r)
    
      print(str_before, file=text_file)
    i+=1
    j+=1


  i = 999
  # j is counting stitches per round divided by 10
  j = 1
  for r in pra[:-1]:
    #print(r)
    print("j : {}, r/sr: {}".format(j, r/sr))
    angles = np.linspace(0, 2*np.pi, j*10+1)
    k = 0
    if (j < int(annotate_index)/2+1):
      print("Printing Red Ring {}".format(j))
      colore="#ff0000"
    else:
      colore="#ffffff"

    for a in angles:
      mod_factor = 1
      #if (k % 6 == 0):
      #  mod_factor = 2
      if (k % j == 0):
        mod_factor = 4
      if (k == 0):
        mod_factor = 8

      # check that the distance from each eye bounding circle center to end point is less than R
      xo = cx+(j-1)*sr*np.cos(a)
      yo = cy+(j-1)*sr*np.sin(a)
      xp = sr*np.cos(a)
      yp = sr*np.sin(a)

      cy = 747.37-(747.37-715.59)
      cxl = 311.26001
      cxr = 2411.26
      # oouter-end
      dl = np.sqrt((xo+xp-cxl)**2 + (yo+yp-cy)**2)
      dr = np.sqrt((xo+xp-cxr)**2 + (yo+yp-cy)**2)
      # inner-end
      dli = np.sqrt((xo-cxl)**2 + (yo-cy)**2)
      dri = np.sqrt((xo-cxr)**2 + (yo-cy)**2)
      
      if (dl > R) | (dr > R):
        #if starting points are interior, then shorten to intersection
        if (dli < R) & (dri < R):
          print("we're interior! {} < {} & {} < {}".format(dli, R, dri, R))
          # are we on the left or right?
          if (a >= 0) & (a < np.pi/2):
            # angle of reference is pi minus angle?
            aa = np.pi - a
          elif (a >= np.pi/2) & (a < np.pi):
            aa = a
          elif (a >= np.pi) & (a < 3*np.pi/2):
            # angle of reference is pi minus angle?
            aa = 2*np.pi - a
          elif (a >= 3*np.pi/2) & (a < 2*np.pi):
            aa = a - np.pi

          theta_l = np.arcsin(L/R*np.sin(aa))
          theta_r = np.pi-aa-theta_l
          print("angles: theta_a: {}, theta_l: {}, theta_r: {}".format(aa, theta_l, theta_r))
          rp = np.sqrt(R**2 + L**2 - 2*L*R*np.cos(theta_r))
          dr = r - rp
          print("sr = {}, dr: {} = {} - {} = r - rp".format(sr, dr, r, rp))
          xo = xo -dr*np.cos(a)#+ (sr-dr)*np.cos(a)
          yo = yo -dr*np.sin(a)#+ (sr-dr)*np.sin(a)
          xp = (dr+sr)*np.cos(a)
          yp = (dr+sr)*np.sin(a)
                  
        str_next = "<path\n     style=\"fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n     d=\"m {},{} {},{} z\"\n     id=\"path{}\"\n     inkscape:connector-curvature=\"0\" />".format(colore, lwidth*mod_factor, xo, yo, xp, yp, i)
        print(str_next, file=text_file)
      i+=1
      k+=1
    j+=1


  print(str_end, file=text_file)
