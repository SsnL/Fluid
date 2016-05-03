#!/usr/bin/python
import sys
import xml.etree.cElementTree as ET

'''
copy and paste from http://effbot.org/zone/element-lib.htm#prettyprint
it basically walks your tree and adds spaces and newlines so the tree is
printed in a nice way
'''
def indent(elem, level=0):
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i

def add_particle(ps, pos, v):
  p = ET.SubElement(ps, "particle")
  ET.SubElement(p, "pos").text = str(pos[0])+" "+str(pos[1])+" "+str(pos[2])
  ET.SubElement(p, "v").text = str(v[0])+" "+str(v[1])+" "+str(v[2])

def build_tree(size, dis, density):
  particles = ET.Element("particles")
  # density = ET.SubElement(particles, "density").text = str(density)#rest density
  density = ET.SubElement(particles, "density").text = str(150.0)#rest density
  ps = ET.SubElement(particles, "ps")

  #fs = ET.SubElement(particles, "fs")

  # for i in range(1-size, size):
  #   for j in range(0,size):
  #     for k in range(1-size,size):
  #       add_particle(ps, [dis*i,dis*j+0.5,dis*k], [0,0,0])


    # for (int i = -6; i < 6; i++)
    #   for (int j = 0; j < 5; j++)
    #     for (int k = -6; k < 6; k++)
    #       ps.push_back(new Particle(
    #         Vector3D(0.15 * i, 0.2 * j + 0.5, 0.15 * k),
    #         Vector3D(0, -0.01, 0),
    #         150.0f
    #       ));
  for i in range(-6, 6):
    for j in range(0, 5):
      for k in range(-6, 6):
        add_particle(ps, [0.15 * i, 0.2 * j + 0.5, 0.15 * k], [0, -0.01, 0])

  indent(particles)
  return ET.ElementTree(particles)

args = sys.argv
if len(args) != 5:
  print "  usage: ./pgen.py <size> <dis> <rest_density> <outputfile>"
  print "  Example: ./pgen.py 3 0.05 1.0 p.xml"
  print "  <size> controls the number of particles"
  print "  <dis> controls how close the particles are to each other"
else:
  size = int(args[1])
  dis  = float(args[2])
  density = float(args[3])
  outputfile = str(args[4])

  tree = build_tree(size,dis,density)
  tree.write(outputfile)
