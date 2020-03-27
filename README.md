# Setup environment
Just to create a conda environment using the Requirements.txt file executing the command as it is explained in the file:
```
# Create the env (naming it 'covid19' as an example)
conda create --name covid19 --file Requirements.txt

# Activate the environment
conda activate covid19

# Run main.py
python src/main.py
``` 

To run the fill pipeline (download sequences, align them, create the tree, visualize it) you have to un-comment that file. 
Keep in min though that aligning the sequences (specially for the 20 countries) will take a couple of hours, depending on your computer.

# Feedback
Feel free to reach out if you have any comment, suggestion or recommendation. Use this code as much as you can and more importantly, in
the way you want.