#' Function to Calculate Cross-Validated R2
#'
#' A function to calculate cross-validated R2
#' @param cv.result
#' @keywords cross-validation R2
#' @export
#' @examples 

PLSK.cv.r2 <- function(cv.result){
  obs <- cv.result$cv.rawdata[,c("native_id","pollutant_conc")]
  names(obs) <- c("native_id","y_obs")
  pred <- cbind.data.frame(names(cv.result$cv.pred),cv.result$cv.pred)
  names(pred) <- c("native_id","y_pred")
  temp <- merge(obs, pred, by="native_id")
  cv.r2(temp$y_obs, temp$y_pred*temp$y_pred)
}

#' Cross-Validation Plot
#'
#' Plot of predictions for out-of-sample locations from cross-validation models
#' against observations.
#' @param cv.result
#' @param ... plotting parameters
#' @keywords cross-validation
#' @export
#' @examples 

cv.plot <- function(cv.result, ...){
  obs <- cv.result$cv.rawdata[,c("native_id","pollutant_conc")]
  names(obs) <- c("native_id","y_obs")
  pred <- cbind.data.frame(names(cv.result$cv.pred),cv.result$cv.pred)
  names(pred) <- c("native_id","y_pred")
  lur.pred    <- cbind.data.frame(names(cv.result$cv.lur),cv.result$cv.lur)
  names(lur.pred) <- c("native_id","lur_pred")
  temp1 <- merge(obs, pred, by="native_id")
  temp  <- merge(temp1, lur.pred, by="native_id")
  temp$y_pred <- temp$y_pred*temp$y_pred
  temp$lur    <- temp$lur*temp$lur
  par(pty="s")
  plot(temp$y_obs, temp$y_pred, 
       xlim=range(c(temp$y_obs, temp$y_pred)), 
       ylim=range(c(temp$y_obs, temp$y_pred)), 
       main=round(summary(lm(temp$y_pred~temp$y_obs))$r.squared, 2),...)
  abline(0,1,col="gray")
  abline(lm(temp$y_pred~temp$y_obs), col="red")
  print(cv.r2(temp$y_obs, temp$y_pred))
  dev.new()
  par(pty="s")
  plot(temp$y_obs, temp$lur, 
       xlim=range(c(temp$y_obs, temp$lur)), 
       ylim=range(c(temp$y_obs, temp$lur)), 
       main=round(summary(lm(temp$lur~temp$y_obs))$r.squared, 2),...)
  abline(0,1,col="gray")
  abline(lm(temp$lur~temp$y_obs), col="red")
  print(cv.r2(temp$y_obs, temp$lur))
}

#' Prediction Plot
#'
#' Plot of predictions against observations.
#' @param full.model
#' @param rawdata
#' @param ... plotting parameters
#' @keywords prediction
#' @export
#' @examples 

pred.plot <- function(full.model, rawdata=full.model$rawdata, ...){
  preds <- PLSK.predict(full.model, rawdata, desc.vars=c("county","state","state_plane","lambert_x","lambert_y"))

  obs  <- full.model$rawdata[,c("location_id","pollutant_conc")]
  names(obs) <- c("location_id","y_obs")
  pred <- preds[,c("location_id", "pred_modelingscale")]
  names(pred) <- c("location_id", "y_pred")
  pred$y_pred <- pred$y_pred*pred$y_pred

  temp <- merge(obs, pred, by="location_id")
  par(pty="s")
  plot(temp$y_obs, temp$y_pred, 
       xlim=range(c(temp$y_obs, temp$y_pred), na.rm=TRUE), 
       ylim=range(c(temp$y_obs, temp$y_pred), na.rm=TRUE), ...)
  abline(0,1,col="gray")
}

#' R2 utility function
#'
#' Returns lm, R2, and RMSE for 2 vectors
#' @param x numeric column
#' @param y numeric column
#' @keywords R2 RMSE
#' @export
#' @examples 

cv.r2 <- function(x, y){
  xy.linear <- lm(y~x)
  xy.r2     <- summary(xy.linear)$r.squared
  xy.rmse   <- sqrt(mean(xy.linear$residuals*xy.linear$residuals))
  list("linear"=xy.linear, "xy.r2"=xy.r2, "xy.rmse"=xy.rmse)
}