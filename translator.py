def trans(arr,vec):
	arr[0][0]+=vec[0]
	arr[1][0]+=vec[0]
	arr[2][0]+=vec[0]
	
	arr[0][1]+=vec[1]
	arr[1][1]+=vec[1]
	arr[2][1]+=vec[1]

	arr[0][2]+=vec[2]
	arr[1][2]+=vec[2]
	arr[2][2]+=vec[2]

def output(arr):
	print "triangle "+str(arr[0][0])+' '+str(arr[0][1])+' '+str(arr[0][1])
	for i in range(1,5):
			print str(arr[i][0])+' '+str(arr[i][1])+' '+str(arr[i][2])
	print  "D"

arr1=[[1,2.745,11.35] ,[2.2,2.745,11.35] ,[1 ,2.745, 10.15] ,[1, 1 ,1] ,[0.5 ,0.5, 0.5]]
arr2=[[2,2.745,11.35] ,[2.2,2.745,11.35] ,[1 ,2.745, 10.15] ,[1, 1 ,1] ,[0.5 ,0.5, 0.5]]
vec=[3,0,0]

# translate 3 in x

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr2)


# translate 3 in x

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in x

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in y
vec=[0,3,0]

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in y

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in y

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate -3 in x
vec=[-3,0,0]

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate -3 in x

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate -3 in x

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in y
vec=[0,-3,0]

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in y

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)

# translate 3 in y

trans(arr1,vec)
output(arr1)

trans(arr2,vec)
output(arr1)
