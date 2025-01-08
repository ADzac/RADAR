from diffimg import diff


file = "Before.jpg"
fileN = "BeforeN.jpg"
fileDB = "BeforeDB.jpg"
file_2 = "New.jpg"
file_2N = "NewN.jpg"
file_2DB = "NewDB.jpg"


'''
print("Figure",round(diff("Realiste.jpg","Realistic.jpg")*100,2))

print("Linear",round(diff("RealisteN.jpg","RealisticN.jpg")*100,2))

print("dB ", round(diff("RealisteDB.jpg","RealisticDB.jpg")*100,2))
'''

#print("Figure",round(diff(file,file_2)*100,2))

print("Linear",round(diff(fileN,file_2N)*100,2))

print("dB ", round(diff(fileDB,file_2DB)*100,2))
