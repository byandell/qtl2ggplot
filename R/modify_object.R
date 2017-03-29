# Various objects have attributes that are worth preserving. This does that.
modify_object <- function(object, new_object) {
  x_attr <- attributes(object)
  x_class <- class(object)
  attrs <- names(x_attr)
  attrs <- attrs[!(attrs %in% c("class", "names", "dim", "dimnames"))]
  
  for(obj in attrs) {
    attr(new_object, obj) <- x_attr[[obj]]
  }
  class(new_object) <- x_class
  new_object
}
