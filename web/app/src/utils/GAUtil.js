import { useEffect } from "react";
import ReactGA from "react-ga";
import { useLocation } from "react-router-dom";

const PROD = (process.env.NODE_ENV ?? "production") === "production";

function GAUtil() {
  const location = useLocation();

  useEffect(() => {
    if (PROD) {
      ReactGA.initialize("UA-162676656-3");
    }
  }, []);

  useEffect(() => {
    if (PROD) {
      ReactGA.pageview(`${location.pathname}${location.search}`);
    }
  }, [location.pathname, location.search]);

  return null;
}

export default GAUtil;
