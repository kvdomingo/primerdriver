import axios from "axios";

const { NODE_ENV } = process.env;

const baseURL = NODE_ENV === "production" ? "https://primerdriver.herokuapp.com/api/" : "http://localhost:8000/api/";

const axiosInstance = axios.create({ baseURL });

const api = {
  data: {
    primerDriver(body) {
      return axiosInstance.post("/", body);
    },
    version() {
      return axiosInstance.get("/version");
    },
    expressionSystem() {
      return axiosInstance.get("/expressionsys");
    },
  },
};

export default api;
