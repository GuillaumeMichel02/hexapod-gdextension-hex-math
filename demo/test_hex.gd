extends Node2D


func _ready():
    var dist = HexMath.distance(Vector2i(0,2), Vector2i(2,1))
    print(dist)  # Should print 3
    var line = Vector2i(-1,4)
    var coord = Vector2i(2,1) # Setup time measurement
    print(HexMath.distance_to_line(coord, line, Vector2i(0,0)))
    
