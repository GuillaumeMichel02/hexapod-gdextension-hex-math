extends Node2D


func _ready():
    var dist = HexMath.distance(Vector2i(-2,0), Vector2i(0,2))
    print(dist)  # Should print 3
    var from = Vector2i(0, 0)
    var to = Vector2i(3920, 6400)
    # Setup time measurement
    var start_time = Time.get_ticks_msec()
    HexMath.get_line(from, to)
    var end_time = Time.get_ticks_msec()
    print("Time taken for get_line: ", end_time - start_time, " ms")
     
    start_time = Time.get_ticks_msec()
    HexMath.hex_dda(from, to)  
    end_time = Time.get_ticks_msec()
    print("hex_dda time taken: ", end_time - start_time, " ms")
