import axios from "axios";

const baseURL = "/api";

const axiosInstance = axios.create({ baseURL });

const api = {
  data: {
    primerDriver(body) {
      return axiosInstance.post("", body);
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
