# Assign value using <- to a name (variable)
x <- 5          # numeric scalar
name <- "ASP"   # character scalar, don't forget the quotes
#this is a comment, do not execude this y<-3


#Vectors: values of the same type
abundance <- c(90, 85, 78)   # numeric vector
names <- c("Bacillus","Ecoli", "Strep") # character vector

#Data frames: table like structure with rows/column
#columns can be different types
df <- data.frame(
  Name = c("Bacillus","Ecoli", "Strep"),
  Abundance = c(90,85,78) 
)

#Accessing elements
x           # scalar
abundance[2]   # 2nd element in vector → 85
df$Abundance    # column from dataframe
df[1,2]     # row 1, column 2 → 90

#Functions: perform tasks on data 
# function_name(arguments)
mean(abundance) #average of numeric vector
sum(abundance)  #sum of scores
