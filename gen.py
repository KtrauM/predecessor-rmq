import random
N = 10
print(N)
for i in range(N):
    # print(random.randint(-50,50))
    print(random.randint(-9223372036854775808, 9223372036854775807))
# for i in range(N):
#     for j in range(i, N):
#         print(f'{i},{j}')

for i in range(N):
    k = random.randint(0, N)
    l = random.randint(0, N)
    if k > l:
        k, l = l, k
    print(f'{k},{l}')